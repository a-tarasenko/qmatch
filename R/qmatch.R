#'Qmatch algorithm
#'
#'Given treatment and control groups, finds a pairing of treatment units
#'to controls that maximizes number of matched disjoint pairs satisfying a
#'given caliper.
#'
#'Qmatch is a generic function, calls one of its methods.
#'
#'@param x Numeric vector with the scores of treated and control objects,
#'           or formula, or glm.
#'@param z If x is a numeric vector, z is a vector with 1 for treated objects
#            and 0 for controls, otherwise z is a dummy argument.
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
#'
#'@return A qmatch object containing:
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
#' @return\code{discarded.t}:   vector of the numbers (positions in \code{x})
#'                of treated objects which are not matched.
#' @return\code{discarded.c}:   vector of the numbers (positions in \code{x})
#'                of control objects which are not matched.
#' @return\code{total.matches}:  the number of matched treated objects, which is the
#'                  number of rows in match.matrix.
#' @return\code{total.pairs}:   the number of matched control objects, which is the
#'                number of non-NA elements in columns \code{2:(controls+1)}
#'                of match.matrix
#'
#'@author Pavel S. Ruzankin, Marina V. Muravleva
#'
#'@references
#' P.S. Ruzankin (2019). A fast algorithm for maximal propensity score matching.
#' \emph{Methodol Comput Appl Probab}. \url{https://doi.org/10.1007/s11009-019-09718-4}
#'
#'
#'@export
#'@docType methods
#'@rdname qmatch
#'@examples qmatch(c(1,1.1,2.5,1.5,0.2,2,0.5,2),c(0,1,0,1,0,1,0,1),0.5)

qmatch <- function(x,
                   z,
                   caliper,
                   controls=1L,
                   data,
                   within,
                   method="nno",
                   m.order="largest")
{
  if (missing(data)) {
    UseMethod("qmatch")
  } else {
    classvar <- vector()
    class(classvar) <- eval(substitute(class(x)),data)
    UseMethod("qmatch",classvar)
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
#   discarded.t: vector of the numbers (positions in scores.t)
#                of treated objects which are not matched.
#   discarded.c: vector of the numbers (positions in scores.c)
#                of control objects which are not matched.
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
  match.matrix <- matrix(as.numeric(NA), nrow = len.scores.t, ncol = 1L+controls)
  discarded.t <- numeric(len.scores.t)
  discarded.c <- numeric(len.scores.c)
  num_pairs <- integer(len.scores.t)
  total.pairs <- 0 # current number of matched controls (of pairs)
  total.matches <- 0 # current number of matched treated


  match.matrix[,1L]=seq_len(len.scores.t)



  # The first element of the right pointers vector is the pointer from outer space
  # Therefore always add 1 to the index!!!
  pointr.c <- c(seq_along(scores.c),0)

  # The last element of the left pointers vector is the pointer from outer space
  pointl.c <- 0:length(scores.c)


  for (control in seq_len(controls)) {

    curc <- pointr.c[1]
    m <- 0 #Number of matches for this control number
    #Vectors for temporary storage of matched pairs for each control nnumber
    match.t <- numeric(len.scores.t)
    match.c <- numeric(len.scores.t)

    for (curt in seq_along(scores.t)) {

      #First, pass as many controls as needed
      while (curc != 0 && scores.c[curc]<scores.t[curt]) {
        curc <- pointr.c[curc+1]
      }

      if (curc != 0) {
        prevc=pointl.c[curc]
        if (prevc != 0 ) {
          #Select the nearest control and check the caliper
          if (scores.t[curt]-scores.c[prevc] >
              scores.c[curc]-scores.t[curt]) {
            if (scores.c[curc]-scores.t[curt] <= caliper) {
              m <- m+1
              match.t[m] <- curt
              match.c[m] <- curc
              nextc <- pointr.c[curc+1]
              pointr.c[prevc+1] <- nextc
              if (nextc != 0) {
                pointl.c[nextc] <- prevc
              } else {
                pointl.c[len.scores.c+1] <- prevc
              }
              curc <- nextc
            }
          } else {
            if (scores.t[curt]-scores.c[prevc] <=  caliper) {
              m <- m+1
              match.t[m] <- curt
              match.c[m] <- prevc
              prevprevc=pointl.c[prevc]
              pointr.c[prevprevc+1] <- curc
              pointl.c[curc] <- prevprevc
            }
          }
        } else { # prevc == 0
          if (scores.c[curc]-scores.t[curt] <= caliper) {
            m <- m+1
            match.t[m] <- curt
            match.c[m] <- curc
            nextc <- pointr.c[curc+1]
            pointr.c[prevc+1] <- nextc
            if (nextc != 0) {
              pointl.c[nextc] <- prevc
            } else {
              pointl.c[len.scores.c+1] <- prevc
            }
            curc <- nextc
          }
        }
      } else { # curc == 0
        prevc <- pointl.c[len.scores.c+1]
        if (prevc != 0) {
          if (scores.t[curt]-scores.c[prevc] <=  caliper) {
            m <- m+1
            match.t[m] <- curt
            match.c[m] <- prevc
            prevprevc=pointl.c[prevc]
            pointr.c[prevprevc+1] <- curc
            pointl.c[len.scores.c+1] <- prevprevc
          }
        }
      }

    }

    if (m==0) break

    total.pairs <- total.pairs + m

    match.t <- head(match.t,m)
    num_pairs[match.t] <- num_pairs[match.t] + 1L

    match.c <- head(match.c,m)
    match.c <- match.c[order(scores.c[match.c])] # Optimal rematching

    match.matrix[match.t,1L+control] <- match.c

  }

  selectrows <- num_pairs > 0L
  discarded.t <- which(!selectrows)

  match.matrix <- matrix(match.matrix[selectrows,],ncol=1L+controls)
  # matrix() used to hadle cases with one or no rows selected

  num_pairs <- num_pairs[selectrows]
  total.matches <- length(num_pairs)

  selcontrols <- rep(TRUE,len.scores.c)
  selcontrols[as.vector(match.matrix[,2L:(1L+controls)])] <- FALSE
  discarded.c <- which(selcontrols)

  list(match.matrix = match.matrix, num_pairs = num_pairs,
       discarded.t = discarded.t, discarded.c = discarded.c,
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
#   discarded.t: vector of the numbers (positions in scores.t)
#                of treated objects which are not matched.
#   discarded.c: vector of the numbers (positions in scores.c)
#                of control objects which are not matched.
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
  match.matrix <- matrix(as.integer(NA), nrow = len.scores.t, ncol = 1L+controls)
  discarded.t <- integer(len.scores.t)
  discarded.c <- integer(len.scores.c)
  num_pairs <- integer(len.scores.t)
  total.pairs <- 0L # current number of matched controls (of pairs)
  total.matches <- 0L # current number of matched treated
  i <- 1L # current control object
  j <- 1L # current treated object
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
      discarded.c[i-total.pairs] <- i
      i <- i + 1L
    } else {
      if (k == 0L) {
        discarded.t[j-total.matches] <- j
      }
      k <- 0L
      j <- j + 1L
    }
  }

  if (k > 0L) {
    j <- j + 1L
  }
  if (j <= len.scores.t) {
    discarded.t[j:len.scores.t - total.matches] <- j:len.scores.t
  }
  if (i <= len.scores.c) {
    discarded.c[i:len.scores.c - total.pairs] <- i:len.scores.c
  }
  match.matrix <- head(match.matrix, n = total.matches)
  discarded.t <- head(discarded.t, n = len.scores.t-total.matches)
  discarded.c <- head(discarded.c, n = len.scores.c-total.pairs)
  num_pairs <- head(num_pairs, n = total.matches)

  list(match.matrix = match.matrix, num_pairs = num_pairs,
       discarded.t = discarded.t, discarded.c = discarded.c,
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
#   discarded.t: vector of the numbers (positions in x)
#                of treated objects which are not matched.
#   discarded.c: vector of the numbers (positions in x)
#                of control objects which are not matched.
#   total.matches: the number of matched treated objects, which is the
#                  number of rows in match.matrix.
#   total.pairs: the number of matched control objects, which is the
#                number of non-NA elements in cols 2:(control+1)
#                of match.matrix
#' @export
qmatch.numeric <- function(x,
                           z,
                           caliper,
                           controls=1L,
                           data,
                           within,
                           method="nno",
                           m.order="largest")
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
      pos.t <- which(z!=0L)
      pos.c <- which(z==0L)
      scores.t <- x[pos.t] # unordered scores for treated
      scores.c <- x[pos.c] # unordered scores for controls
    } else {
      pos.t <- eval(substitute(which(z!=0L)),data)
      pos.c <- eval(substitute(which(z==0L)),data)
      scores.t <- eval(substitute(x[z!=0L]),data) # unordered scores for treated
      scores.c <- eval(substitute(x[z==0L]),data) # unordered scores for controls
    }
    
    if (m.order == "largest") {
      scores.t <- -scores.t
      scores.c <- -scores.c
    }
    perm.t <- order(scores.t)
    perm.c <- order(scores.c)

    scores.t <- scores.t[perm.t] # order scores for treated
    scores.c <- scores.c[perm.c] # order scores for controls
    if (tolower(method)=="qmatch") {
      qm <- qmatch_core(scores.t, scores.c, caliper, controls)
    } else if (tolower(method)=="nno") {
      qm <- nnomatch_core(scores.t, scores.c, caliper, controls)
    } else {
      stop(paste0("Wrong method ",method))
    }

    # Replace positions in scores.t and scores.c with positions in x
    qm$match.matrix[,1L] <- pos.t[perm.t[qm$match.matrix[,1L]]]
    qm$match.matrix[,2L:(controls+1L)] <-
      pos.c[perm.c[qm$match.matrix[,2L:(controls+1L)]]]
    qm$discarded.t <- pos.t[perm.t[qm$discarded.t]]
    qm$discarded.c <- pos.c[perm.c[qm$discarded.c]]

    class(qm) <- "qmatch"
    qm
  } else { #if (missing(within))
    if (missing(data)){
      if (length(within)!=length(z)) {
        stop("lengths of within and z differ")
      }

      fact <- unique(within)
      match.matrix <- matrix(as.integer(NA), nrow =length(x),
          ncol = 1L+controls)
      discarded.t <- integer(length(x))
      discarded.c <- integer(length(x))
      num_pairs <- integer(length(x))
      total.pairs <- 0L
      total.matches <- 0L
      total.discarded.t <- 0L
      total.discarded.c <- 0L

      for (m in seq_along(fact)) {
        currselect <- which(within==fact[m])
        qm <- qmatch.numeric(x[currselect], z[currselect],
            caliper, controls, method=method, m.order=m.order)

        if (qm$total.matches>0L) {
          match.matrix[total.matches + 1L:qm$total.matches,] <-
              currselect[qm$match.matrix]
          num_pairs[total.matches + 1L:qm$total.matches] <- qm$num_pairs
        }
        total.matches <- total.matches + qm$total.matches
        total.pairs <- total.pairs +  qm$total.pairs
        if (length(qm$discarded.c)>0L) {
          discarded.c[total.discarded.c + seq_along(qm$discarded.c)] <-
              currselect[qm$discarded.c]
          total.discarded.c <- total.discarded.c + length(qm$discarded.c)
        }
        if (length(qm$discarded.t)>0L) {
          discarded.t[total.discarded.t + seq_along(qm$discarded.t)] <-
              currselect[qm$discarded.t]
          total.discarded.t <- total.discarded.t + length(qm$discarded.t)
        }
      }
      match.matrix <- head(match.matrix, n = total.matches)
      discarded.t <- head(discarded.t, n = total.discarded.t)
      discarded.c <- head(discarded.c, n = total.discarded.c)
      num_pairs <- head(num_pairs, n = total.matches)

    } else { # if (missing(data))

      if (eval(substitute(length(within)),data)!=eval(substitute(length(z)),data)) {
        stop("lengths of within and z differ")
      }

      fact <- unique(eval(substitute(within),data))
      match.matrix <- matrix(as.integer(NA),
          nrow =eval(substitute(length(x)),data), ncol = 1L+controls)
      discarded.t <- integer(eval(substitute(length(x)),data))
      discarded.c <- integer(eval(substitute(length(x)),data))
      num_pairs <- integer(eval(substitute(length(x)),data))
      total.pairs <- 0L
      total.matches <- 0L
      total.discarded.t <- 0L
      total.discarded.c <- 0L
      for (m in seq_along(fact)) {
        pos.t <- eval(substitute(which(z!=0L&(within==fact[m]))),data)
        pos.c <- eval(substitute(which(z==0L&(within==fact[m]))),data)
        scores.t <- eval(substitute(x[z!=0L&(within==fact[m])]),data) # unordered scores for treated
        scores.c <- eval(substitute(x[z==0L&(within==fact[m])]),data) # unordered scores for controls

        if (m.order == "largest") {
          scores.t <- -scores.t
          scores.c <- -scores.c
        }
        perm.t <- order(scores.t)
        perm.c <- order(scores.c)

        scores.t <- scores.t[perm.t] # order scores for treated
        scores.c <- scores.c[perm.c] # order scores for controls

        if (tolower(method)=="qmatch") {
          qm <- qmatch_core(scores.t, scores.c, caliper, controls)
        } else if (tolower(method)=="nno") {
          qm <- nnomatch_core(scores.t, scores.c, caliper, controls)
        } else {
          stop(paste0("Wrong method ",method))
        }

        #Replace positions in scores.t and scores.c with positions in x
        qm$match.matrix[,1L] <- pos.t[perm.t[qm$match.matrix[,1L]]]
        qm$match.matrix[,2L:(controls+1L)] <-
          pos.c[perm.c[qm$match.matrix[,2L:(controls+1L)]]]
        qm$discarded.t <- pos.t[perm.t[qm$discarded.t]]
        qm$discarded.c <- pos.c[perm.c[qm$discarded.c]]


        if (qm$total.matches>0L) {
          match.matrix[total.matches + 1L:qm$total.matches,] <- qm$match.matrix
          num_pairs[total.matches + 1L:qm$total.matches] <- qm$num_pairs
        }
        total.matches <- total.matches + qm$total.matches
        total.pairs <- total.pairs +  qm$total.pairs
        if (length(qm$discarded.t)>0L) {
          discarded.t[total.discarded.t + seq_along(qm$discarded.t)] <- qm$discarded.t
          total.discarded.t <- total.discarded.t + length(qm$discarded.t)
        }
        if (length(qm$discarded.c)>0L) {
          discarded.c[total.discarded.c + seq_along(qm$discarded.c)] <- qm$discarded.c
          total.discarded.c <- total.discarded.c + length(qm$discarded.c)
        }
      }
      match.matrix <- head(match.matrix, n = total.matches)
      discarded.t <- head(discarded.t, n = total.discarded.t)
      discarded.c <- head(discarded.c, n = total.discarded.c)
      num_pairs <- head(num_pairs, n = total.matches)
    }
    qm <- list(match.matrix = match.matrix, num_pairs = num_pairs,
                 discarded.t = discarded.t, discarded.c = discarded.c,
                 total.matches = total.matches, total.pairs = total.pairs)
    class(qm) <- "qmatch"
    qm

  }
}

# qmatch.formula: decomposes the formula
#                     and calls qmatch.numeric function
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
#   discarded.t: vector of the numbers (positions in x)
#                of treated objects which are not matched.
#   discarded.c: vector of the numbers (positions in x)
#                of control objects which are not matched.
#   total.matches: the number of matched treated objects, which is the
#                  number of rows in match.matrix.
#   total.pairs: the number of matched control objects, which is the
#                number of non-NA elements in cols 2:(control+1)
#
#' @export
qmatch.formula <- function(x,
                           z,
                           caliper,
                           controls=1L,
                           data,
                           within,
                           method="nno",
                           m.order="largest")
{
  if (length(x)!=3L || x[[2]]==".") {
    stop("Formula must have a left hand side")
  }
  if (x[[3]]==".") {
    stop("Formula must have a right hand side")
  }
  if (length(x[[2]])>2L) {
    stop("Left part of the formula is complex")
  }
  if (length(x[[3]])>2L) {
    stop("Right part of the formula is complex")
  }
  # Below unclass() is needed when I() function is used in the formula.
  if (missing(within)) {
    if (missing(data)) {
      ax <- eval(unclass(x[[3]]), envir = parent.frame())
      az <- eval(unclass(x[[2]]), envir = parent.frame())
    } else {
      ax <- eval(unclass(x[[3]]), data)
      az <- eval(unclass(x[[2]]), data)
    }
    qmatch.numeric(ax,az,caliper,controls,method=method,m.order=m.order)
  } else {
    if (missing(data)) {
      ax <- eval(unclass(x[[3]]), envir = parent.frame())
      az <- eval(unclass(x[[2]]), envir = parent.frame())
      awithin <- within
    } else {
      ax <- eval(unclass(x[[3]]), data)
      az <- eval(unclass(x[[2]]), data)
      awithin <- eval(substitute(within),data)
    }
    qmatch.numeric(ax,az,caliper,controls,within=awithin,
                   method=method,m.order=m.order)
  }
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
#   discarded.t: vector of the numbers (positions in x)
#                of treated objects which are not matched.
#   discarded.c: vector of the numbers (positions in x)
#                of control objects which are not matched.
#   total.matches: the number of matched treated objects, which is the
#                  number of rows in match.matrix.
#   total.pairs: the number of matched control objects, which is the
#                number of non-NA elements in cols 2:(control+1)
#
#' @export
qmatch.glm <- function(x,z,
                       caliper,
                       controls=1L,
                       data,
                       within,
                       method="nno",
                       m.order="largest")
{
  if (!missing(data)) {
    stop("The data argument cannot be used to specify the environment for glm")
  }
  if (length(x$na.action) >= 1L) {
    warning("There were missing values in the glm model. The resulting match.matrix will contain references to vectors with excluded missing values")
  }
  if (missing(within)) {
    qmatch.numeric(x$fitted.values,x$model[,1],caliper,controls,
                   method=method,m.order=m.order)
  } else {
    qmatch.numeric(x$fitted.values,x$model[,1],caliper,controls,
                   within = eval(substitute(within),x$data), 
                   method=method,m.order=m.order)
  }
}

#' @export
qmatch.default <- function(x,...)
{
  stop(paste0("Class \"",class(x),"\" is not supported for qmatch()"))
}


