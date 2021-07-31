#'Two algorithms for quick matching without replacement on a scalar index.
#'The first algorithm is greedy nearest neighbor matching 
#'followed by optimal rematching.
#'The second algorithm finds a pairing of treatment units
#'to controls that maximizes the number of matched pairs satisfying the
#'given caliper.
#'
#'The input data may be presented as numeric vector \code{x} of scores and 
#'vector \code{z} representing the comparison group. 
#'Observations where \code{z!=0} are assigned 
#'to the treatment group, the other observations being assigned to the control 
#'group. Besides, the \code{x} argument can be specified as a formula 
#'\code{group~score}. Or \code{x} may be a \code{glm} object, 
#'\code{fitted.values} of which are treated as scores and 
#'the dependent variable is treated as group.
#'
#'
#'@param x Numeric vector with the scores of treated and control objects,
#'           or formula, or glm.
#'@param z If \code{x} is a numeric vector, then \code{z} is a vector with nonzero 
#'           vales for treated objects
#            and zeros for controls, otherwise \code{z} is a dummy argument. 
#'@param caliper The caliper, i.e., the maximal distance allowed between
#            the scores of matched treated and control object.
#'@param controls The maximal number of controls matched to a single treated object
#            (the minimal number is always 1 for the algorithm).
#'@param data An optional data frame, list or environment for arguments x,z
#'           and within.
#'@param within Vector of factors dividing objects into groups.
#'@param method Method for matching: "nno" (default) for nearest neighbor matching 
#'           followed by optimal rematching, or "qmatch" for the method producing 
#'           maximal number of pairs under the given caliper.
#'@param m.order Matching order: "largest" for matching in descending order 
#'           starting from largest values of the score, 
#'           "smallest" for matching in ascending order starting from 
#'           smallest values of the score.
#'@param compute.weights Logical value indicating whether to compute weights. 
#'@param compute.cluster Logical value indicating whether to compute cluster IDs.
#'@param compute.discarded Logical value indicating whether to compute 
#'           numbers of discarded observations.
#'
#'@return A \code{qmatch} object containing:
#' @return\code{match.matrix}:   a \code{(controls+1)} columns matrix, each row containing the numbers
#'                 of a matched treated (in col 1) and corresponding
#'                 control (in columnss \code{2:(controls+1))} objects. The numbers
#'                 of objects are positions in vector x.
#'                 If, for a treated object, there are only \code{0<k<controls}
#'                 matches, the elements of the row in the
#'                 cols (k+2):(controls+1) are filled with integer NAs.
#' @return\code{num_pairs}:   vector of length total.matches showing the number of
#'              controls matched to each treated object. \code{num_pairs[j]}
#'              corresponds to the j-th row of match.matrix.
#' @return\code{total.matches}:  the number of matched treated objects, which is the
#'                  number of rows in match.matrix.
#' @return\code{total.pairs}:   the number of matched control objects, which is the
#'                number of non-NA elements in columns \code{2:(controls+1)}
#'                of match.matrix
#' @return\code{weights}:   Vector of weights of the observations. 
#'                If an observation is not matched, then the corresponding element of 
#'                \code{weights} is 0. If a treated observation is matched, 
#'                then the corresponding element of \code{weights} is 1.
#'                If a control observation is matched, then the corresponding element of
#'                \code{weights} is 1 over the number of controls matched 
#'                to the corresponding treated observation.
#' @return\code{cluster}:   Vector of cluster IDs. Each cluster consists of 
#'                a treated observation and all the control observations matched to
#'                the treated one. 
#'                If an observation is not matched, then the corresponding element of 
#'                \code{cluster} is NA.
#' @return\code{discarded}:   Vector of the numbers (positions in \code{x})
#'                of observations that are not matched.
#'
#'
#'@author Pavel S. Ruzankin, Marina V. Muravleva
#'
#'@references
#' P.S. Ruzankin (2020) A fast algorithm for maximal propensity score matching.
#' \emph{Methodol Comput Appl Probab}. \url{https://doi.org/10.1007/s11009-019-09718-4}
#'
#'
#'@export
#'@docType methods
#'@rdname qmatch
#'@examples qmatch(c(1,1.1,2.5,1.5,0.2,2,0.5,2),c(0,1,0,1,0,1,0,1),0.5)

qmatch <- function(x, z, caliper, controls = 1, data, within, method = "nno", m.order = "largest", compute.weights = TRUE, compute.cluster = TRUE, compute.discarded = TRUE)
{
  if (missing(data)) {
    UseMethod("qmatch")
  } else {
    classvar <- vector()
    class(classvar) <- eval(substitute(class(x)), data)
    UseMethod("qmatch", classvar)
  }
}


# nnomatch_core: core function for optimized NNM algorithm
#
# Arguments:
#   scores.t: sorted vector with the scores of treated objects.
#   scores.c: sorted vector with the scores of control objects.
#   caliper:  the caliper, i.e., the maximal distance allowed between
#             the scores of matched treated and control object.
#   controls: the maximal number of controls matched to a single treated object
#             (the minimal number is always 1 for the algorithm).
#
#
# Value: list with the following elements.
#   match.matrix: a (control+1) col matrix, each row containing the numbers
#                 of a matched treated (in col 1) and corresponding
#                 control (in cols 2:(control+1)) objects. The numbers
#                 of objects are positions in vectors scores.t and scores.c.
#                 If, for a treated object, there are only 0<k<controls
#                 matches, the elements of the row in the
#                 cols (k+2):(controls+1) are filled with integer NAs.
#   num_pairs: vector of length total.matches showing the number of
#              controls matched to each treated object. num_pairs[j]
#              corresponds to the j-th row of match.matrix.
#   total.matches: the number of matched treated objects, which is the
#                  number of rows in match.matrix.
#   total.pairs: the number of matched control objects, which is the
#                number of non-NA elements in cols 2:(control+1)
#                of match.matrix


nnomatch_core <- function(scores.t,
                          scores.c,
                          caliper,
                          controls = 1L)
{
  controls <- as.integer(controls)
  len.scores.c <- length(scores.c)
  len.scores.t <- length(scores.t)
  type.scores <- typeof(max(len.scores.c,len.scores.t))
  # as.vector() is usually more efficient than as(); it does not matter here
  match.matrix <- matrix(as.vector(NA, mode = type.scores), 
                         nrow = len.scores.t, ncol = 1L + controls)
  num_pairs <- integer(len.scores.t)
  total.pairs <- as.vector(0L, mode = type.scores)  # current number of matched controls (of pairs)
  total.matches <- as.vector(0L, mode = type.scores) # current number of matched treated

  match.matrix[,1L] <- seq_len(len.scores.t)

  # The first element of the right pointers vector is the pointer from outer space
  # Therefore always add 1 to the index!!!
  pointr.c <- c(seq_along(scores.c), 0L)

  # The last element of the left pointers vector is the pointer from outer space
  pointl.c <- 0L:length(scores.c)


  for (control in seq_len(controls)) {

    cur.c <- pointr.c[1L]
    m <- as.vector(0L, mode = type.scores) # Number of matches (matched treated) for this control number iteration
    # Vectors for temporary storage of matched pairs for each control number
    match.t <- vector(mode = type.scores, len.scores.t)
    matched.c <- rep.int(FALSE, len.scores.c)

    for (cur.t in seq_along(scores.t)) {

      #First, pass as many controls as needed
      while (cur.c != 0L && scores.c[cur.c]<scores.t[cur.t]) {
        cur.c <- pointr.c[cur.c + 1L]
      }

      if (cur.c != 0L) {
        prev.c=pointl.c[cur.c]
        if (prev.c != 0L ) {
          #Select the nearest control and check the caliper
          if (scores.t[cur.t] - scores.c[prev.c] >
              scores.c[cur.c] - scores.t[cur.t]) {
            if (scores.c[cur.c] - scores.t[cur.t] <= caliper) {
              m <- m + 1L
              match.t[m] <- cur.t
              matched.c[cur.c] <- TRUE
              next.c <- pointr.c[cur.c + 1L]
              pointr.c[prev.c + 1L] <- next.c
              if (next.c != 0L) {
                pointl.c[next.c] <- prev.c
              } else {
                pointl.c[len.scores.c + 1L] <- prev.c
              }
              cur.c <- next.c
            }
          } else {
            if (scores.t[cur.t] - scores.c[prev.c] <=  caliper) {
              m <- m + 1L
              match.t[m] <- cur.t
              matched.c[prev.c] <- TRUE
              prevprev.c <- pointl.c[prev.c]
              pointr.c[prevprev.c + 1L] <- cur.c
              pointl.c[cur.c] <- prevprev.c
            }
          }
        } else { # prev.c == 0
          if (scores.c[cur.c] - scores.t[cur.t] <= caliper) {
            m <- m + 1L
            match.t[m] <- cur.t
            matched.c[cur.c] <- TRUE
            next.c <- pointr.c[cur.c + 1L]
            pointr.c[prev.c + 1L] <- next.c
            if (next.c != 0L) {
              pointl.c[next.c] <- prev.c
            } else {
              pointl.c[len.scores.c + 1L] <- prev.c
            }
            cur.c <- next.c
          }
        }
      } else { # cur.c == 0
        prev.c <- pointl.c[len.scores.c + 1L]
        if (prev.c != 0L) {
          if (scores.t[cur.t] - scores.c[prev.c] <=  caliper) {
            m <- m + 1L
            match.t[m] <- cur.t
            matched.c[prev.c] <- TRUE
            prevprev.c <- pointl.c[prev.c]
            pointr.c[prevprev.c + 1L] <- cur.c
            pointl.c[len.scores.c + 1L] <- prevprev.c
          }
        }
      }

    }

    if (m == 0L) break

    total.pairs <- total.pairs + m

    match.t <- head(match.t, m)
    num_pairs[match.t] <- num_pairs[match.t] + 1L

    match.matrix[match.t, 1L + control] <- which(matched.c) # Optimal rematching

  }

  selectrows <- num_pairs > 0L

  match.matrix <- matrix(match.matrix[selectrows,], ncol = 1L + controls)
  # here matrix() is used to handle cases with one or no rows selected

  num_pairs <- num_pairs[selectrows]
  total.matches <- length(num_pairs)

  list(match.matrix = match.matrix, num_pairs = num_pairs,
       total.matches = total.matches, total.pairs = total.pairs)
}



# qmatch_core: core function for qmatch algorithm
#
# Arguments:
#   scores.t: sorted vector with the scores of treated objects.
#   scores.c: sorted vector with the scores of control objects.
#   caliper:  the caliper, i.e., the maximal distance allowed between
#             the scores of matched treated and control object.
#   controls: the maximal number of controls matched to a single treated object
#             (the minimal number is always 1 for the algorithm).
#
# Value: list with the following elements.
#   match.matrix: a (control+1) col matrix, each row containing the numbers
#                 of a matched treated (in col 1) and corresponding
#                 control (in cols 2:(control+1)) objects. The numbers
#                 of objects are positions in vectors scores.t and scores.c.
#                 If, for a treated object, there are only 0<k<controls
#                 matches, the elements of the row in the
#                 cols (k+2):(controls+1) are filled with integer NAs.
#   num_pairs: vector of length total.matches showing the number of
#              controls matched to each treated object. num_pairs[j]
#              corresponds to the j-th row of match.matrix.
#   total.matches: the number of matched treated objects, which is the
#                  number of rows in match.matrix.
#   total.pairs: the number of matched control objects, which is the
#                number of non-NA elements in cols 2:(control+1)
#                of match.matrix

qmatch_core <- function(scores.t,
                        scores.c,
                        caliper,
                        controls = 1L)
{
  controls <- as.integer(controls)
  len.scores.c <- length(scores.c)
  len.scores.t <- length(scores.t)
  type.scores <- typeof(max(len.scores.c,len.scores.t))
  match.matrix <- matrix(as.vector(NA, mode = type.scores), 
                         nrow = len.scores.t, ncol = 1L + controls)
  num_pairs <- integer(len.scores.t)
  total.pairs <- as.vector(0L, mode = type.scores) # current number of matched controls (of pairs)
  total.matches <- as.vector(0L, mode = type.scores) # current number of matched treated
  i <- as.vector(1L, mode = type.scores) # current control object
  j <- as.vector(1L, mode = type.scores) # current treated object
  k <- 0L # current number of controls matched to the current treated
  while (i <= len.scores.c && j <= len.scores.t) {
    if (abs(scores.c[i] - scores.t[j]) <= caliper) {
      total.pairs <- total.pairs + 1L
      k <- k + 1L
      if (k == 1L) {
        total.matches <- total.matches + 1L
        match.matrix[total.matches, 1L] <- j
      }
      num_pairs[total.matches] <- k
      match.matrix[total.matches, k+1L] <- i
      i <- i + 1L
      if (k >= controls) {
        k <- 0L
        j <- j + 1L
      }
    } else if (scores.c[i] < scores.t[j]) {
      i <- i + 1L
    } else {
      k <- 0L
      j <- j + 1L
    }
  }

  if (k > 0L) {
    j <- j + 1L
  }
  match.matrix <- head(match.matrix, n = total.matches)
  num_pairs <- head(num_pairs, n = total.matches)

  list(match.matrix = match.matrix, num_pairs = num_pairs,
       total.matches = total.matches, total.pairs = total.pairs)
}

# qmatch.numeric: sorts the vector x and calls qmatch_core function
#
# Arguments:
#   x: vector with the scores of treated and control objects
#   z: vector with 1 for treated objects and 0 for controls.
#   caliper: the caliper, i.e., the maximal distance allowed between
#            the scores of matched treated and control object.
#   data: an optional data frame, list or environment containing the
#         vectors x,z and within.
#   within: vector of factors for exact matching (stratification).
#   method: method for matching: "nno" for optimized NNM or "qmatch".
#   m.order: order of matching: begin from "largest" or from "smallest".
#   compute.weights: whether to compute weights. 
#   compute.cluster: whether to compute cluster IDs.
#   compute.discarded: whether to compute numbers of discarded observations.
#
# Value: object of class qmatch containing:
#   match.matrix: a (control+1) col matrix, each row containing the numbers
#                 of a matched treated (in col 1) and corresponding
#                 control (in cols 2:(control+1)) objects. The numbers
#                 of objects are positions in vector x.
#                 If, for a treated object, there are only 0<k<controls
#                 matches, the elements of the row in the
#                 cols (k+2):(controls+1) are filled with integer NAs.
#   num_pairs: vector of length total.matches showing the number of
#              controls matched to each treated object. num_pairs[j]
#              corresponds to the j-th row of match.matrix.
#   total.matches: the number of matched treated objects, which is the
#                  number of rows in match.matrix.
#   total.pairs: the number of matched control objects, which is the
#                number of non-NA elements in cols 2:(control+1)
#                of match.matrix
#   discarded:   vector of the numbers (positions in x) of observations
#                that are not matched.

#' @export
qmatch.numeric <- function(x,
                           z,
                           caliper,
                           controls=1L,
                           data,
                           within,
                           method="nno",
                           m.order="largest",
                           compute.weights = TRUE,
                           compute.cluster = TRUE,
                           compute.discarded = TRUE)
{
  controls <- as.integer(controls)

  if (controls < 1L) {
    stop("controls < 1")
  }
  if (missing(caliper)){
    stop("caliper is not declared")
  }
  
  m.order <- tolower(m.order)
  if ( ! (m.order %in% c("largest","smallest")) ) {
    stop('m.order can only be "largest" or "smallest"')
  }

  if (missing(data)){
    if (anyNA(x)) {
      stop("NA in x")
    }
    if (anyNA(z)) {
      stop("NA in z")
    }
    if (length(x)!=length(z)) {
      stop("lengths of x and z differ")
    }
  } else {
    if (eval(substitute(anyNA(x)),data)) {
      stop("NA in x")
    }
    if (eval(substitute(anyNA(z)),data)) {
      stop("NA in z")
    }
    if (eval(substitute(length(x)),data)!=eval(substitute(length(z)),data)) {
      stop("lengths of x and z differ")
    }
  }

  if (missing(within)) {
    if (missing(data)) {
      len.scores <- length(x) # needed for weights calculation
      
      pos.t <- which(z!=0L)
      pos.c <- which(z==0L)
      scores.t <- x[pos.t] # unordered scores for treated
      scores.c <- x[pos.c] # unordered scores for controls
    } else {
      len.scores <- eval(substitute(length(x)), data) # needed for weights calculation
      
      pos.t <- eval(substitute(which(z!=0L)), data)
      pos.c <- eval(substitute(which(z==0L)), data)
      scores.t <- eval(substitute(x[z!=0L]), data) # unordered scores for treated
      scores.c <- eval(substitute(x[z==0L]), data) # unordered scores for controls
    }
    
    if (m.order == "largest") {
      scores.t <- -scores.t
      scores.c <- -scores.c
    }
    perm.t <- order(scores.t)
    perm.c <- order(scores.c)

    scores.t <- scores.t[perm.t] # order scores for treated
    scores.c <- scores.c[perm.c] # order scores for controls
    if (tolower(method) == "qmatch") {
      qm <- qmatch_core(scores.t, scores.c, caliper, controls)
    } else if (tolower(method) == "nno") {
      qm <- nnomatch_core(scores.t, scores.c, caliper, controls)
    } else {
      stop(paste0("Wrong method ",method))
    }

    # Replace positions in scores.t and scores.c with positions in x
    qm$match.matrix[,1L] <- pos.t[perm.t[qm$match.matrix[,1L]]]
    qm$match.matrix[,2L:(controls + 1L)] <-
      pos.c[perm.c[qm$match.matrix[,2L:(controls + 1L)]]]

  } else { #if (missing(within))
    
    if (missing(data)){
      if (length(within) != length(x)) {
        stop("lengths of within and x differ")
      }

      len.scores <- length(x) # needed for weights calculation
      
      fact <- unique(within)
      type.scores <- typeof(length(x))
      match.matrix <- matrix(as.vector(NA, mode = type.scores), 
                             nrow =length(x),
                             ncol = 1L+controls)
      num_pairs <- integer(length(x))
      total.pairs <- as.vector(0L, mode = type.scores)
      total.matches <- as.vector(0L, mode = type.scores)

      for (m in seq_along(fact)) {
        currselect <- which(within == fact[m])
        qm <- qmatch.numeric(x[currselect], z[currselect],
            caliper, controls, method = method, m.order = m.order,
            compute.weights = FALSE,
            compute.cluster = FALSE,
            compute.discarded = FALSE)

        if (qm$total.matches > 0L) {
          match.matrix[total.matches + 1L:qm$total.matches,] <-
              currselect[qm$match.matrix]
          num_pairs[total.matches + 1L:qm$total.matches] <- qm$num_pairs
        }
        total.matches <- total.matches + qm$total.matches
        total.pairs <- total.pairs +  qm$total.pairs
      }
      match.matrix <- head(match.matrix, n = total.matches)
      num_pairs <- head(num_pairs, n = total.matches)

    } else { # if (missing(data))

      if (eval(substitute(length(within)),data) != eval(substitute(length(z)),data)) {
        stop("lengths of within and z differ")
      }

      fact <- unique(eval(substitute(within), data))
      len.scores <- eval(substitute(length(x)), data)
      type.scores <- typeof(len.scores)
      match.matrix <- matrix(as.vector(NA, mode = type.scores),
          nrow = len.scores, ncol = 1L + controls)
      num_pairs <- integer(len.scores)
      total.pairs <- as.vector(0L, mode = type.scores)
      total.matches <- as.vector(0L, mode = type.scores)
      for (m in seq_along(fact)) {
        pos.t <- eval(substitute(which(z!=0L&(within==fact[m]))), data)
        pos.c <- eval(substitute(which(z==0L&(within==fact[m]))), data)
        scores.t <- eval(substitute(x[z!=0L&(within==fact[m])]), data) # unordered scores for treated
        scores.c <- eval(substitute(x[z==0L&(within==fact[m])]), data) # unordered scores for controls

        if (m.order == "largest") {
          scores.t <- -scores.t
          scores.c <- -scores.c
        }
        perm.t <- order(scores.t)
        perm.c <- order(scores.c)

        scores.t <- scores.t[perm.t] # order scores for treated
        scores.c <- scores.c[perm.c] # order scores for controls

        if (tolower(method) == "qmatch") {
          qm <- qmatch_core(scores.t, scores.c, caliper, controls)
        } else if (tolower(method) == "nno") {
          qm <- nnomatch_core(scores.t, scores.c, caliper, controls)
        } else {
          stop(paste0("Wrong method ",method))
        }

        #Replace positions in scores.t and scores.c with positions in x
        qm$match.matrix[,1L] <- pos.t[perm.t[qm$match.matrix[,1L]]]
        qm$match.matrix[,2L:(controls + 1L)] <-
          pos.c[perm.c[qm$match.matrix[,2L:(controls + 1L)]]]


        if (qm$total.matches > 0L) {
          match.matrix[(total.matches + 1L):(total.matches + qm$total.matches),] <- 
            qm$match.matrix
          num_pairs[(total.matches + 1L):(total.matches + qm$total.matches)] <- 
            qm$num_pairs
        }
        total.matches <- total.matches + qm$total.matches
        total.pairs <- total.pairs +  qm$total.pairs
      }
      match.matrix <- head(match.matrix, n = total.matches)
      num_pairs <- head(num_pairs, n = total.matches)
    } # if (missing(data)) else
    
    qm <- list(match.matrix = match.matrix, num_pairs = num_pairs,
               total.matches = total.matches, total.pairs = total.pairs)
  } # if (missing(within)) else
  
  
  if (compute.weights) {
    if (controls == 1L) {
      weights <- rep.int(0L, len.scores)
      weights[qm$match.matrix] <- 1L
    } else {
      weights <- rep.int(0, len.scores)
      weights[qm$match.matrix[,1L]] <- 1
      
      for (j in seq_len(qm$total.matches)) {
        weights[qm$match.matrix[j, 2L:(1L + controls)]] <- 1 / (qm$num_pairs[j])
      }
    }
    qm$weights <- weights
  }
  
  if (compute.cluster) {
    cluster <- rep.int(as.vector(NA, mode = typeof(nrow(qm$match.matrix))), 
                       len.scores)
    for (j in seq_len(qm$total.matches)) {
      cluster[qm$match.matrix[j,]] <- j
    }
    qm$cluster <- cluster
  }
  
  if (compute.discarded) {
    select <- rep.int(TRUE, len.scores)
    select[qm$match.matrix] <- FALSE
    qm$discarded <- which(select)
  }
  
  
  class(qm) <- "qmatch"
  qm
  
}

# qmatch.formula: decomposes the formula
#                 and calls qmatch.numeric function
#
# Arguments:
#   x: formula of the form az~ax, where ax is the vector with the scores
#      of treated and control objects, and az is the treatment indicator
#      (the vector with 1 for treated objects and 0 for controls).
#   z: dummy argument.
#   caliper: the caliper, i.e., the maximal distance allowed between
#            the scores of matched treated and control object.
#   controls: the maximal number of controls matched to a single treated object
#             (the minimal number is always 1 for the algorithm).
#   data: an optional data frame, list or environment containing the
#         vectors ax and az.
#   method: method for matching: "nno" for optimized NNM or "qmatch".
#   m.order: order of matching: begin from "largest" or from "smallest".
#   compute.weights: whether to compute weights. 
#   compute.cluster: whether to compute cluster IDs.
#   compute.discarded: whether to compute numbers of discarded observations.
#
# Value: object of class qmatch containing:
#   match.matrix: a (control+1) col matrix, each row containing the numbers
#                 of a matched treated (in col 1) and corresponding
#                 control (in cols 2:(control+1)) objects. The numbers
#                 of objects are positions in vector ax.
#                 If, for a treated object, there are only 0<k<controls
#                 matches, the elements of the row in the
#                 cols (k+2):(controls+1) are filled with integer NAs.
#   num_pairs: vector of length total.matches showing the number of
#              controls matched to each treated object. num_pairs[j]
#              corresponds to the j-th row of match.matrix.
#   total.matches: the number of matched treated objects, which is the
#                  number of rows in match.matrix.
#   total.pairs: the number of matched control objects, which is the
#                number of non-NA elements in cols 2:(control+1)
#   discarded:   vector of the numbers (positions in x) of observations
#                that are not matched.
#
#' @export
qmatch.formula <- function(x,
                           z,
                           caliper,
                           controls = 1L,
                           data,
                           within,
                           method = "nno",
                           m.order = "largest",
                           compute.weights = TRUE,
                           compute.cluster = TRUE,
                           compute.discarded = TRUE)
{
  if (length(x) != 3L || x[[2L]] == ".") {
    stop("Formula must have a left hand side")
  }
  if (x[[3L]] == ".") {
    stop("Formula must have a right hand side")
  }
  if (length(x[[2L]]) > 2L) {
    stop("Left part of the formula is complex")
  }
  if (length(x[[3L]]) > 2L) {
    stop("Right part of the formula is complex")
  }
  # Below unclass() is needed when I() function is used in the formula.
  if (missing(within)) {
    if (missing(data)) {
      ax <- eval(unclass(x[[3L]]), envir = parent.frame())
      az <- eval(unclass(x[[2L]]), envir = parent.frame())
    } else {
      ax <- eval(unclass(x[[3L]]), data)
      az <- eval(unclass(x[[2L]]), data)
    }
    qmatch.numeric(ax, az, caliper, controls, 
                   method = method, m.order = m.order,
                   compute.weights = compute.weights,
                   compute.cluster = compute.cluster,
                   compute.discarded = compute.discarded)
  } else { # if (missing(within))
    if (missing(data)) {
      ax <- eval(unclass(x[[3L]]), envir = parent.frame())
      az <- eval(unclass(x[[2L]]), envir = parent.frame())
      awithin <- within
    } else {
      ax <- eval(unclass(x[[3L]]), data)
      az <- eval(unclass(x[[2L]]), data)
      awithin <- eval(substitute(within), data)
    }
    qmatch.numeric(ax, az, caliper, controls, within=awithin,
                   method = method, m.order = m.order,
                   compute.weights = compute.weights,
                   compute.cluster = compute.cluster,
                   compute.discarded = compute.discarded)
  } # if (missing(within)) else
}


# qmatch.glm: decomposes the glm variable
#                 and calls qmatch.numeric function
#
# Arguments:
#   x: glm varibale with the formula of the form az~ax,
#      where ax is the vector with the scores of treated and control objects,
#      and az is the treatment indicator
#      (the vector with 1 for treated objects and 0 for controls),
#      and an optional data frame, list or environment for the formula.
#   z: dummy argument.
#   caliper: the caliper, i.e., the maximal distance allowed between
#            the scores of matched treated and control object.
#   controls: the maximal number of controls matched to a single treated object
#             (the minimal number is always 1 for the algorithm).
#   data: dummy argument.
#   method: method for matching: "nno" for optimized NNM or "qmatch".
#   m.order: order of matching: begin from "largest" or from "smallest".
#   compute.weights: whether to compute weights. 
#   compute.cluster: whether to compute cluster IDs.
#   compute.discarded: whether to compute numbers of discarded observations.
#
# Value: object of class qmatch containing:
#   match.matrix: a (control+1) col matrix, each row containing the numbers
#                 of a matched treated (in col 1) and corresponding
#                 control (in cols 2:(control+1)) objects. The numbers
#                 of objects are positions in vector x.
#                 If, for a treated object, there are only 0<k<controls
#                 matches, the elements of the row in the
#                 cols (k+2):(controls+1) are filled with integer NAs.
#   num_pairs: vector of length total.matches showing the number of
#              controls matched to each treated object. num_pairs[j]
#              corresponds to the j-th row of match.matrix.
#   total.matches: the number of matched treated objects, which is the
#                  number of rows in match.matrix.
#   total.pairs: the number of matched control objects, which is the
#                number of non-NA elements in cols 2:(control+1)
#   discarded:   vector of the numbers (positions in x) of observations
#                that are not matched.
#
#' @export
qmatch.glm <- function(x,z,
                       caliper,
                       controls = 1L,
                       data,
                       within,
                       method = "nno",
                       m.order = "largest",
                       compute.weights = TRUE,
                       compute.cluster = TRUE, 
                       compute.discarded = TRUE)
{
  if (!missing(data)) {
    stop("The data argument cannot be used to specify the environment for glm")
  }
  if (length(x$na.action) >= 1L) {
    warning("There were missing values in the glm model. The resulting match.matrix will contain references to vectors with excluded missing values")
  }
  if (missing(within)) {
    qmatch.numeric(x$fitted.values, x$model[,1L], caliper, controls,
                   method = method, m.order = m.order,
                   compute.weights = compute.weights,
                   compute.cluster = compute.cluster,
                   compute.discarded = compute.discarded)
  } else {
    qmatch.numeric(x$fitted.values, x$model[,1L], caliper, controls,
                   within = eval(substitute(within),x$data), 
                   method = method, m.order = m.order,
                   compute.weights = compute.weights,
                   compute.cluster = compute.cluster,
                   compute.discarded = compute.discarded)
  }
}

#' @export
qmatch.default <- function(x, ...)
{
  stop(paste0("Class \"", class(x), "\" for x is not supported in qmatch()"))
}
