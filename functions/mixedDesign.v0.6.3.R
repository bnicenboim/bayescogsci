########  GENERATE A DATAFRAME WITH SIMULATED DATA ###########################################
# Sven Hohenstein, Reinhold Kliegl, 2009-2015

# Version 0.6.3 (June 2015)
#
# changes:
# - removed all `cat`s

################
# INPUT:
#
# B:      vector of between-subject factors levels (it's length equals the number of between-subject factors)
#         (you can use B = c() or B = NULL if there are no between-subject factors)
# W:      vector of within-subject factors levels (it's length equals the number of within-subject factors)
#         (you can use W = c() or W = NULL if there are no within-subject factors)
# n:      number of observations within each between-subjects factor level combination (natural number)
#         (the total number of observations equals n * prod(B); the minimum n equals prod(W), but can not be lower than 2)
# M:	  Matrix specifying the cell means of crossing between- and within-subject factors (dimensions: prod(B)-by-prod(W) matrix)
#         (for pure within-subjects designs, it is possible to input a vector of size prod(W) as an alternative to a 1-by-prod(W) matrix)
#         OR
#         a scalar (single number) that will be the mean of all cells
# SD:     Matrix specifying the cell standard deviations of crossing between- and within-subject factors (dimensions: prod(B)-by-prod(W) matrix)
#         (for pure within-subjects designs, it is possible to input a vector of size prod(W) as an alternative to a 1-by-prod(W) matrix)
#         OR
#         a scalar (single number) that will be the standard deviation of all cells
# R:	  List of matrices specifying the correlations of the witin-subject factor level combinations
#		  (list length: prod(B); dimensions of matrices: prod(W)-by-prod(W))
#         OR
#		  a singe matrix that will be used for the correlation of all within-subject factor level combinations   
#		  OR      
#		  a scalar (single number) that will be the correlation of all within-subject factor level combinations (except the diagonal)
#		  OR 
#		  a list of length prod(B) containing only single numbers or single numbers and matrices
#         Note: R is is ignored if there are no within-subject factors
# miss:   Matrix specifying the proportion of random data loss within each cell of crossing between- and within-subject factors
#		  (dimensions: prod(B)-by-prod(W) matrix)
#         (for pure within-subjects designs, it is possible to input a vector of size prod(W) as an alternative to a 1-by-prod(W) matrix)
#         OR
#         a scalar (single number), [0, 1], that will be the relative amount of data loss in each cell
# long:   if TRUE each line of the resulting data frame contains of one value of the dependent variable
#         if FALSE all (within-subject) measures of one subject are in one line of the resulting data frame
# family: string specifying the output distribution; "gaussian", "gamma", or "beta"
# empirical: if TRUE the random values will match the specified properties exactly (at least for gaussian data)
#            if FALSE all values will be drawn randomly from a distribution with the specified properties
################

################
# OUTPUT:
#
# data frame containing the labels for the between- and-within factors, the subject (id) numbers, and the measurements obtained randomly from a normal distribution
# (the values of the dependent variable are random BUT they match the mean, standard deviation, and correlation values of the input)
################

################
# Please note:
# This works perfectly for guassian distributions.
# But note two suboptimal instances occuring with gamma and beta distributions:
# 	1) Obtained means and standard deviations will slightly differ from the input. The deviation decreases with the number of subjects (observations).
#	2) Obtained correlations will slightly differ from the input. The deviation decreases with the number of subjects (observations).
#	   Whereas the classical Pearson correlation statistic is informative for gaussian distributions, it is hardly convenient for gamma and beta distributions.
#	   If a gamma or beta distribution is used, the Spearman rank correlation is the appropriate statistic. This method is used if cor() is called with the parameter method="spearman".
#	   The obtained Spearman correlations will be more close to the input values than Pearson correlations.
################


mixedDesign <- function(B = NULL, W = NULL, n = ifelse(is.null(W), 2, prod(W)), M = 0, SD = 1, R = 0, miss = 0, family = "gaussian", long = FALSE, empirical = TRUE) {
  if (class(try(library(MASS), silent = FALSE)) == "try-error") stop("MASS package is necessary for this function. Download and install it.")
  ##### catch wrong input, e.g., B <- c(1, 2, 3)
  if (any(B < 2)) {	# remove illegal conditions (less than 2 conditions per factor)
    warning("Removing illegal between-subject factors.")
    ix <- which(B < 2)
    B <- B[- ix] 	
  }
  nB <- length(B)
  if (any(W < 2)) {	# remove illegal conditions (less than 2 conditions per factor)
    warning("Removing illegal within-subject factors.")
    ix <- which(W < 2)
    W <- W[- ix] 	
  }
  nW <- length(W) # number of within-subject factors
  if (nB || nW) {
    nBcomb <- 1	# number of within-subject factor level combinations
    # this start value remains constant (1) for cases of pure within-subject designs since there is no "X by 0" design
    if (nB) {	# this segment is ignored if there are no between-subject factors
      Blabel <- NULL # list with item labels
      for (i in 1:nB) {
      	Blabel[[i]] <- paste(LETTERS[i], 1:B[i], sep = "")
      	nBcomb <- nBcomb * B[i]
      }
    }
    nWcomb <- 1	# number of within-subject factor level combinations
    # this start value remains constant (1) for cases of pure between-subject designs since there is no "X by 0" design
    if (nW) {	# this segment is ignored if there are no within-subject factors
      Wlabel <- NULL # list with item labels
      Wlabel_comb <- NULL # combination of factor level items
      for (i in 1:nW) {
        Wlabel_comb <- rep(Wlabel_comb, each = W[[i]])  
        Wlabel[[i]] <- paste(letters[i], 1:W[i], sep = "")
        if (i == 1) Wlabel_comb <- Wlabel[[i]] else Wlabel_comb <- paste(Wlabel_comb, rep(Wlabel[[i]], nWcomb), sep = "_")  # arrangement of combinations
        nWcomb <- nWcomb * W[i]
      }
      # add "W" to the labels of the within factors
      Wlabel_comb <- paste("W_", Wlabel_comb, sep ="") 
    }
    if (n < 2) stop("Not enough subjects for the given design")
    if (n < nWcomb) {	# if subjects are not sufficient to produce data frame, additional subjects are added and later removed
    	warning("Not enough subjects for the given design:\n the resulting means, standard deviations, and correlations may differ from the specified ones")
    	n.dif <- nWcomb - n	# subjects to remove
    	n <- nWcomb    	
   	} else n.dif <- 0	# no subjects to remove
	# Means: allow vectors as an alternative to 1-by-X matrices
    if (!nB && !is.matrix(M)) {
      M <- as.matrix(M) # necessary for comparing dimensions
      dim(M) <- rev(dim(M)) # change vector to 1-by-X matrix 
    }
    if (!all(dim(as.matrix(M)) == c(nBcomb, nWcomb))) {
      #generate matrix of mean values)
      if (!all(dim(as.matrix(M)) == c(1, 1))) warning("Inanpropriate matrix of means. Using first value for all cell means.") # more than one element
      CellMeans <- matrix(M[1], nBcomb, nWcomb)
    } else CellMeans <- M
    # SDs: allow vectors as an alternative to 1-by-X matrices
    if (!nB && !is.matrix(SD)) {
      SD <- as.matrix(SD) # necessary for comparing dimensions
      dim(SD) <- rev(dim(SD)) # change vector to 1-by-X matrix 
    }
    if (!all(dim(as.matrix(SD)) == c(nBcomb, nWcomb))) {
      #generate matrix of sd values)
      if (!all(dim(as.matrix(SD)) == c(1, 1))) warning("Inanpropriate dimensions of matrix of SDs. Using first value for all cell SDs.")
      CellSDs <- matrix(SD[1], nBcomb, nWcomb)
    } else CellSDs <- SD 
    if (any(CellSDs < 0)) stop("Negative values are not allowed for standard deviation parameter.")     
    if (!nB && !is.matrix(miss)) {
      miss <- as.matrix(miss) # necessary for comparing dimensions
      dim(miss) <- rev(dim(miss)) # change vector to 1-by-X matrix 
    }
    if (!all(dim(as.matrix(miss)) == c(nBcomb, nWcomb))) {
      #generate matrix of data loss proportions)
      if (!all(dim(as.matrix(miss)) == c(1, 1))) warning("Inanpropriate dimensions of matrix of data loss proportions. Using first value for all cell data loss proportions.")
      Cell.loss <- matrix(miss[1], nBcomb, nWcomb)
    } else Cell.loss <- miss
    # checking data loss values, must be between 0 (no data loss) and 1 (complete data loss)
    if (any(Cell.loss < 0 | Cell.loss > 1)) {
    	stop("Inappropriate values for data loss proportions. All values must be between 0 and 1.")    	}  
    # check correlation matrices
    if (any(unlist(R) < -1 | unlist(R) > 1) & nW) {
    	stop("Inappropriate values for correlation matrix. All values must be between -1 and 1.")    
    } 	 
    N <- n*nBcomb         # total number of subjects
    # ... generate dataframe compatible with means, standard deviation, and cors for w-s factor
    dat.w <- matrix(rnorm(nWcomb * N), N)
    # check correlation matrix
	if (nW) {      
   	  # put single matrix in list (if there are no between-subject factors)    
      if (!is.list(R)) {
    	if (!nB) {
    		R <- list(R)    		     		
    	} else {
    		# use single matrix for all correlation matrices
    		R <- replicate(nBcomb, list(R))		
    		warning("Using identical correlation matrix for all groups.")
    		oneCorMat <- TRUE	# one correlation matrix (or one single number) for all correlations
    	}	
      } else { # test length of correlation list R
      	if (nB & length(R) != nBcomb) stop("Inappropriate matrix list 'R' specified. 'R' must be a list with a length according to the number of (between-subject) groups or a single matrix or a single number.")
      	oneCorMat <- FALSE	# multiple correlation matrices (or single values) 
      }    
      for (group in 1 : nBcomb) {      
      	if (!all(dim(as.matrix(R[[group]])) == c(nWcomb, nWcomb))) {
        	# does the inanpropriate matrix contain more than a single value?
        	if (!all(dim(as.matrix(R[[group]])) == c(1, 1))) {
    			warning(paste("Inanpropriate dimensions of ", group, ". correlation matrix. Using first value for all correlations (except the diagonal).", sep = ""))
        	} else if (nB) if (!oneCorMat) warning("Using ", R[[group]][1], " as value for all correlations in group ", group, ".\n", sep ="") # not if there is a single matrix/value      			
          # create correlation matrix        
        	R[[group]] <- matrix(R[[group]][1], nWcomb, nWcomb)
        	diag(R[[group]]) <- 1
      	} else {
        	# since only the values below the diagonal are relevant, the values above are replaced by them
        	R[[group]][upper.tri(R[[group]])] <- R[[group]][lower.tri(R[[group]])]
        	diag(R[[group]]) <- 1
      	}
      }	  
    } else R <- replicate(nBcomb, list(1)) # if there are no correlations (no within-subject factors)       
    for (group in 1 : nBcomb) {       
      # group index
      ix <- (n * (group - 1) + 1) : (n * group)	
      # function mvrnorm() from MASS package
      dat.w[ix, ] <- mvrnorm(n = n, mu = rep(0, nWcomb), Sigma = R[[group]], empirical = empirical)
    }   
    # remove additional subjects (if n < prod(W))
    if (n.dif != 0) {
    	ix.rem <- NULL	#lines to be removed
    	for (i in 1:nBcomb) {
    		ix <- (n * (i - 1) + 1) : (n * i)	# subjects within one between-subject factor level combination (lines)
	    	ix <- ix[(n - n.dif + 1) : n]	# lines with additional subjects (to be removed)
	    	ix.rem <- c(ix.rem, ix)
    	}
    	dat.w <- dat.w[setdiff(1:N, ix.rem), ] # remove additional subject lines
    	n <- n - n.dif	# give n its actual value back
    	N <- n * prod(B)	# the actual value of N
   	}  
   	# ... between-group effect
    if (nB) {
      #dat.b <- expand.grid(id=1:n, group=Blabel)
      dat.b_small <- expand.grid(rev(Blabel))
      dat.b <- data.frame(matrix(NA, N, nB))
      colnames(dat.b) <- paste("B", LETTERS[1:nB], sep ="_")
        for (i in nB:1) {
          dat.b[ , nB + 1 - i] <- rep(dat.b_small[, i], each = n)
        }
      dat.b$id <- factor(1:N)
      #dat.b$id <- factor(paste(dat.b$group, dat.b$id, sep=""))
    } else {
      dat.b <- data.frame(id = factor(1:N)) # a data frame containing the subject id
    }
    # do operations according distribution family
    if (tolower(family) == "gaussian") {
    	# ... put mean and standard deviation in data
    	for (i in 1:nBcomb) {
      		ix <- (n * (i - 1) + 1) : (n * i)	# subjects within one between-subject factor level combination (lines)
      		for (j in 1:n) {	
        		dat.w[ix[j], ] <- dat.w[ix[j], ] * CellSDs[i, ]
        		dat.w[ix[j], ] <- dat.w[ix[j], ] + CellMeans[i, ]
      		}
    	}    
    } else if (tolower(family) == "gamma") {
    	if (any(CellMeans <= 0)) stop("Inappropriate mean for gamma distribution specified. Choose M > 0.")
    	if (any(CellSDs <= 0)) stop("Inappropriate standard deviation for gamma distribution specified. Choose SD > 0.")
      # generate random beta distributed values
    	for (i in 1:nBcomb) {
      		ix <- (n * (i - 1) + 1) : (n * i)	# subjects within one between-subject factor level combination (lines)
      		for (k in 1:nWcomb)	{
        		# compute rank of gaussian distributed data of one factor level combination (between and within)
        		rg <- rank(dat.w[ix, k])
        		# calculate alpha (shape parameter) and theta (scale parameter) to obtain the desired mean and standard deviation
        		alpha <- CellMeans[i, k]^2 / CellSDs[i, k]^2
        		theta <- CellSDs[i, k]^2 / CellMeans[i, k]
        		# create random gamma data
        		dat.w[ix, k] <- rgamma(n, shape = alpha, scale = theta) 
        		# apply rank of gaussian distributed values to ordered gamma distributed values
				dat.w[ix, k] <- dat.w[ix, k][order(dat.w[ix, k])][rg]
      		}
    	}     
    } else if (tolower(family) == "beta") {
    	if (any(0 >= CellMeans | CellMeans >= 1)) stop("Mean of beta distribution must be in the closed interval (0, 1).")
    	if (any(0 >= CellSDs | CellSDs >= 0.5)) stop("Standard deviation of beta distribution must be in the closed interval (0, 0.5).")
      # generate random beta distributed values
    	for (i in 1:nBcomb) {
      		ix <- (n * (i - 1) + 1) : (n * i)	# subjects within one between-subject factor level combination (lines)
      		for (k in 1:nWcomb)	{
        		# compute rank of gaussian distributed data of one factor level combination (between and within)
        		rg <- rank(dat.w[ix, k])
        		# calculate alpha1 (shape parameter) and alpha2 (shape parameter) to obtain the desired mean and standard deviation
        		alpha <- -((1 / CellSDs[i, k]^2) * CellMeans[i, k] * (CellMeans[i, k]^2 - CellMeans[i, k] + CellSDs[i, k]^2))
        		beta <- (1 / CellSDs[i, k]^2) * (CellMeans[i, k] - 1) * (CellMeans[i, k]^2 - CellMeans[i, k] + CellSDs[i, k]^2)
        		# create random beta data
        		dat.w[ix, k] <- rbeta(n, alpha, beta)
        		if (any(is.na(dat.w[ix, k]))) {
        			stop(paste("Failed to create beta distribution with mean ", CellMeans[i, k], " and standard deviation ", CellSDs[i, k], ". Choose different parameters.", sep = ""))	    		}
        		# apply rank of gaussian distributed values to ordered beta distributed values
				dat.w[ix, k] <- dat.w[ix, k][order(dat.w[ix, k])][rg]
      		}
    	}     
    } else stop("Unknown distribution family; possible are gaussian, gamma, and beta distributions.")   
    # apply data loss randomly
    # 1) change proportion to absolute value
    Cell.loss <- round(Cell.loss * n)
    # 2) apply within each cell of n subjects
    for (i in 1:nBcomb) {
    	ix <- (n * (i - 1) + 1) : (n * i)
    	for (j in 1:nWcomb) {
    		remove <- sample(1:n, Cell.loss[i, j], replace = FALSE)
    		dat.w[ix, j][remove] <- NA
    	}
    }    
    # ... combine dat.b and dat.w  (or dat.l)
    if (long | !nW) {	# there is no wide format for between-subjects designs since each subject is associated with exactly 1 value (DV)
      dat.l <- as.vector(as.matrix(dat.w)) #long format
      dat <- NULL	# this will be the final data frame
      for (i in 1:nWcomb) {
        dat <- rbind(dat, dat.b)
      }
      if (nW) {
        # ... combine with information of within-subject factor levels combinations
        Wcond <- rev(expand.grid(rev(Wlabel)))
        names(Wcond) <- paste("W_", letters[1:nW], sep = "")
        Wcond.l <- NULL
        for (i in 1:N) {  # span a data frame
           Wcond.l <- rbind(Wcond.l, Wcond)
        }
        for (i in 1:nWcomb) {  # modify data frame (each row of Wcond is repeated N times)
          ix <- (N * (i - 1) + 1) : (N * i)
          Wcond.l[ix, ] <- Wcond[i, ]
        }
        dat <- cbind(dat, Wcond.l)	# combine between-subject information and within-subject labels
      }
      dat <- cbind(dat, dat.l)	# combine (between-subject information and within-subject labels) and (measures (DV))
      names(dat)[2 + nB + nW] <- "DV"
      # order data frame along id
      dat <- data.frame(dat[order(dat$id), ], row.names = NULL)
    } else {
      colnames(dat.w) <- Wlabel_comb
      dat <- cbind(dat.b, dat.w)
    }
  } else dat <- NULL	# if there are neither between- nor within-subject factors
return(dat)
} # end of function
