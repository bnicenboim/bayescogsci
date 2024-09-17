# Daniel J. Schad, Maximilian Rabe, Reinhold Kliegl, Feb. 10, 2020

#==============================================================================
simLMM <- function(form, data=NULL, Fixef, RE_sd, RE_cor=0, LMM=NULL, empirical=FALSE, verbose=TRUE, family="gaussian") {
#------------------------------------------------------------------------
# simulate data from a linear mixed-effects model
#------------------------------------------------------------------------

	if (!is.null(LMM)) {

	  if (!class(LMM)%in%c("lmerModLmerTest","lmerMod")) stop('LMM must be a model object from a fit with the lmer() function.')
		EF <- extractLMM(LMM)
		Fixef_ <- EF[[1]]
		Ranef  <- EF[[2]]
		if (is.null(data)) data <- EF[[3]]

	} else {

	  detect_errors1(form, data, RE_sd, RE_cor, family)

	  # extract terms from formula
		Terms <- extractTerms(form, data)

		# collect fixed effects terms (+ fixed-effects estimates)
		Fixef_ <- list(slopes=Terms[[1]], b=Fixef)

		# obtain terms in random effects
		Ranef <- list()
		for (i in 1:length(Terms[[2]])) # i <- 1
			Ranef[[i]] <- list(slopes=Terms[[2]][[i]])

		# expand random effects correlations
		if (!is.list(RE_cor)) {
			RE_cor_list <- list()
			for (i in 1:length(Ranef)) { # i <- 1
				mm_ranef <- model.matrix(as.formula(paste0("~ ",paste0(Terms[[2]][[i]],collapse=" + "))), data=data )
				nre <- dim(mm_ranef)[2]
				#if (Terms[[2]][[i]][1]=="1") nre <- nre - 1
				RE_cor_list[[i]] <- diag(rep(1,nre))
				if (length(RE_cor)==length(Ranef))
					RE_cor_list[[i]][RE_cor_list[[i]]==0] <- RE_cor[i]
				if (length(RE_cor)==1)
					RE_cor_list[[i]][RE_cor_list[[i]]==0] <- RE_cor
			}
			RE_cor <- RE_cor_list
		}

		# Store by, RE sd, and RE correlations
		nvc <- length(Terms[[2]])
		for (i in 1:nvc) { # i <- 1
			Ranef[[i]]$by  <- Terms[[3]][i]
			Ranef[[i]]$sd  <- RE_sd[[i]]
			Ranef[[i]]$cor <- RE_cor[[i]]
		}
		if (family=="gaussian") Ranef[[nvc+1]] <- RE_sd[[nvc+1]]

	}

  detect_errors2(data, Fixef_, Ranef, family, empirical)
  if (empirical==FALSE) resp <- simLMM_empF(data, Fixef_, Ranef, family)
	if (empirical==TRUE ) resp <- simLMM_empT(data, Fixef_, Ranef, family)
  if (verbose  ==TRUE )      simLMM_verbose(data, Fixef_, Ranef, empirical, family)

	return(resp)
}


#==============================================================================
simLMM.b <- function(data, Fixef, Ranef, LMM=NULL, empirical=FALSE) {
#------------------------------------------------------------------------
# simulate data from a linear mixed-effects model
#------------------------------------------------------------------------

	if (!is.null(LMM)) {
		EF <- extractLMM(LMM)
		Fixef <- EF[[1]]
		Ranef <- EF[[2]]
	}

	if (empirical==FALSE) resp <- simLMM_empF(data, Fixef, Ranef)
	if (empirical==TRUE ) resp <- simLMM_empT(data, Fixef, Ranef)

	return(resp)
}


#==============================================================================
detect_errors1 <- function(form, data, RE_sd, RE_cor, family) {

  termlab <- attr(terms(as.formula(form)),"term.labels")
  interc <- attr(terms(as.formula(form)),"intercept")
  # determine random effects
  idx <- c(); for (i in 1:length(termlab)) # i <- 1
    idx[i] <- any(strsplit(termlab[i],"")[[1]] == "|")
  # collect fixed effects terms
  Fixterm <- termlab[!idx]
  if (interc==1) Fixterm <- c("1", Fixterm)
  if (interc==0) Fixterm <- c("0", Fixterm)
  # obtain terms in random effects
  ranterms <- which(idx)
  if (is.list(RE_cor) & length(RE_cor)!=length(ranterms)) stop('Random effects correlations do not fit formula.')
  if (!is.list(RE_cor) & !length(RE_cor)%in%c(1,length(ranterms))) stop('Random effects correlations do not fit formula.')
  if (!is.list(RE_sd)) stop('RE_sd must be specified as a list.')
  if (family=="gaussian") if (length(RE_sd)!=(length(ranterms)+1)) stop('Length of RE_sd does not correspond to formula.')
  if (family!="gaussian") if (length(RE_sd)!= length(ranterms)   ) stop('Length of RE_sd does not correspond to formula.')
  #------ expand formula ------
  for (i in 1:length(ranterms)) { # i <- 1
    formRE <- termlab[ranterms[i]]
    tmp <- which(strsplit(formRE,"")[[1]] == "|")
    formRE1 <- substr(formRE, 1, tmp[1]-1)
    REterms <- attr(terms(as.formula(paste0("~",formRE1))),"term.labels")
    if (length(REterms)!=0) {
      idx3 <- 1:length(REterms)
      if (length(grep(":",REterms))>0) idx3 <- idx3[-grep(":",REterms)]
      if (any(REterms=="0")) idx3 <- idx3[-which(REterms=="0")]
      if (any(REterms=="1")) idx3 <- idx3[-which(REterms=="1")]
      REterms1 <- REterms[idx3]
      if (!all(REterms1 %in% names(data))) stop('Random effects predictor missing in data.')
    }
  }

}

#==============================================================================
detect_errors2 <- function(data, Fixef_, Ranef, family, empirical) {

  data <- data.frame(data)

  # family argument
  if (!any(family%in%c("gaussian","binomial","lp"))) stop('The family argument must be one of "gaussian", "binomial", or "lp" (linear predictor).')
  if( empirical==TRUE & family=="binomial" ) warning('The family argument only partially works with empirical=TRUE.')

  # check whether variables exist in data frame
  # Fixed effects
  idx3 <- 1:length(Fixef_$slopes)
  if (length(grep(":",Fixef_$slopes))>0) idx3 <- idx3[-grep(":",Fixef_$slopes)]
  if (any(Fixef_$slopes=="0")) idx3 <- idx3[-which(Fixef_$slopes=="0")]
  if (any(Fixef_$slopes=="1")) idx3 <- idx3[-which(Fixef_$slopes=="1")]
  fixVars <- Fixef_$slopes[idx3]
  if (!all(fixVars %in% names(data))) stop('Fixed effects predictor missing in data.')
  # Random effects
  if (family=="gaussian") nRanef <- length(Ranef)-1
  if (family!="gaussian") nRanef <- length(Ranef)
  for (i in 1:nRanef) { # i <- 1
    if (!Ranef[[i]]$by %in% names(data)) stop('Grouping variable missing in data.')
    idx3 <- 1:length(Ranef[[i]]$slopes)
    if (length(grep(":",Ranef[[i]]$slopes))>0) idx3 <- idx3[-grep(":",Ranef[[i]]$slopes)]
    if (any(Ranef[[i]]$slopes=="0")) idx3 <- idx3[-which(Ranef[[i]]$slopes=="0")]
    if (any(Ranef[[i]]$slopes=="1")) idx3 <- idx3[-which(Ranef[[i]]$slopes=="1")]
    ranVars <- Ranef[[i]]$slopes[idx3]
    if (!all(ranVars %in% names(data))) stop('Random effects predictor missing in data.')
  }

  # check whether number of coefficients is correct
  # Fixed effects
  mm_fixef <- model.matrix(as.formula(paste0("~ ", paste0(Fixef_$slopes,collapse=" + "))), data=data)
  if (ncol(mm_fixef)!=length(Fixef_$b)) stop('The number of fixed effects coefficients does not match the formula.')
  # Random effects
  if (family!="gaussian") ranef <- Ranef
  if (family=="gaussian") {
    ranef <- list()
    for (i in 1:(length(Ranef)-1)) ranef[[i]] <- Ranef[[i]]
  }
  for (i in 1:length(ranef)) { # i <- 1
    mm_ranef <- model.matrix(as.formula(paste0("~ ",paste0(ranef[[i]]$slopes,collapse=" + "))), data=data )
    if (ncol(mm_ranef)!=length(ranef[[i]]$sd)) stop('The number of random effects coefficients (SDs) does not match the formula.')
    if (ncol(mm_ranef)!=nrow(ranef[[i]]$cor)) stop('The random effects correlation matrix does not match the formula.')
    if (ncol(mm_ranef)!=ncol(ranef[[i]]$cor)) stop('The random effects correlation matrix does not match the formula.')
    if (nrow(ranef[[i]]$cor)!=nrow(ranef[[i]]$cor)) stop('The random effects correlation matrix must have the same number of columns as rows.')
    as.numeric(ranef[[i]]$cor)
  }
}

#==============================================================================
simLMM_verbose <- function(dat, Fixef_, Ranef, empirical, family) {

  if (family=="gaussian")
    cat("Data simulation from a linear mixed-effects model (LMM)\n")
  if (family=="binomial")
    cat("Data simulation from a logistic generalized linear mixed-effects model (GLMM)\n")
  if (family=="lp")
    cat("Data simulation from a linear mixed-effects model (LMM): linear predictor\n")

  Sigma_res <- Ranef[[length(Ranef)]]
  if (family=="gaussian") ranef <- Ranef[1:(length(Ranef)-1)]
  if (family!="gaussian") ranef <- Ranef

  dat <- data.frame(dat)

  # fixed effects: obtain model matrix + write predictor names with coefficients
  if (length(Fixef_$slopes)==1 & Fixef_$slopes[1]=="1") {
    mm_fixef <- model.matrix(as.formula("~ 1"), data=dat)
  } else {
    mm_fixef <- model.matrix(as.formula(paste0("~ ", paste0(Fixef_$slopes,collapse=" + "))), data=dat)
  }
  fixNames <- colnames(mm_fixef)
  bs <- Fixef_$b
  names(bs) <- fixNames
  if (fixNames[1]=="(Intercept)") fixNames[1] <- "1"
  fix_form <- paste0(fixNames,collapse=" + ")
  if (fixNames[1]!="1") fix_form <- paste0("0 + ",fix_form,collapse="")
  #bs  <- as.character(Fixef_$b, length=form_length)


  # random effects: obtain model matrix + compute linear predictor
  form_ranef <- ran_form <- list()
  NC <- ranNames <- ranNames2 <- SDs <- list()
  NCby <- c()
  for (i in 1:length(ranef)) { # i <- 1
    if (length(ranef[[i]]$slopes)==1 & ranef[[i]]$slopes[1]=="1") {
      mm_ranef <- model.matrix(as.formula(paste0("~ 1")), data=dat)
    } else {
      mm_ranef <- model.matrix(as.formula(paste0("~ ",paste0(paste0(ranef[[i]]$slopes,collapse=" + ")))), data=dat )
    }
    ranNames[[i]]  <- colnames(mm_ranef)
    SDs[[i]]       <- ranef[[i]]$sd
    NC[[i]]        <- nchar(ranNames[[i]])
    NCby[i]        <- nchar(ranef[[i]]$by)
    ranNames2[[i]] <- ranNames[[i]]
    if (ranNames2[[i]][1]=="(Intercept)") ranNames2[[i]][1] <- "1"
    if (ranNames2[[i]][1]!="1") ranNames2[[i]][1] <- paste0("0 + ",ranNames2[[i]][1],collapse="")
    ran_form[[i]]  <- paste0("( ", paste0(ranNames2[[i]],collapse=" + "), " | ", ranef[[i]]$by," )")
  }

  FORM <- paste0(fix_form, " + ", paste0(ran_form, collapse=" + "), collapse="")
  cat(paste0("Formula: ",family," ~ "))
  cat(FORM)
  cat("\n")

  nx1 <- max(c(10,max(NCby))) # number of characters for by
  nx2 <- max(c(7,max(sapply(NC,max)))) # number of characters for predictor variables

  #x <- c(100.0001, 1e+1, 1e+0, 1e-1, 1e-2, 1e-3, 1e-10)
  x <- unlist(SDs)
  # number of digits after the comma
  n0d <- floor( -log10( 1e-10 + abs(x) - floor( abs( x ) ) ) )
  n0d[floor(x)!=0] <- 0
  N0D <- min(c(7,max(n0d)+1))
  # number of digits before the comma
  N1D <- max(nchar(as.character(floor(x))))

  if (empirical==TRUE ) cat("empirical = TRUE\n")
  if (empirical==FALSE) cat("empirical = FALSE\n")

  cat("Random effects:\n")
  cat(paste0(" Groups", paste0(rep(" ",nx1-6),collapse=""),
             "Name", paste0(rep(" ",nx2-3),collapse=""),
             "  Std.Dev.", paste0(rep(" ",N0D+N1D-2),collapse=""),
             "Corr\n"))

  for (i in 1:length(ranNames)) { # i <- 1
    for (ii in 1:length(ranNames[[i]])) { # ii <- 1
      # print by name
      if (ii==1) {
        tmp1 <- paste0(" ",ranef[[i]]$by, paste0(rep(" ",nx1-nchar(ranef[[i]]$by)-1),collapse="") )
        nchar(tmp1)
      } else {
        tmp1 <- paste0(rep(" ",nx1),collapse="")
        nchar(tmp1)
      }
      # print predictor variable
      tmp2 <- paste0(" ",ranNames[[i]][ii], paste0(rep(" ",nx2-nchar(ranNames[[i]][ii])),collapse="") )
      # print SD
      tmp3 <- sprintf(paste0("%",N1D+N0D+1,".",N0D,"f"), SDs[[i]][ii])
      # print correlation
      if (ii==1) {
        tmp4 <- ""
      } else {
        #tmp4 <- paste0(sprintf(paste0("%",N1D,".",N0D,"f"), ranef[[i]]$cor[ii,1:(ii-1)]), collapse="   ")
        tmp4 <- paste0(sprintf(paste0("%",N1D,".",2,"f"), ranef[[i]]$cor[ii,1:(ii-1)]), collapse="   ")
      }
      cat(paste0(tmp1,tmp2,"   ",tmp3,"     ",tmp4,"\n",collapse=""))
    }
  }
  if (family=="gaussian")
    cat(paste0(" Residual", paste0(rep(" ",nx1+nx2-5),collapse=""), sprintf(paste0("%",N1D+N0D+1,".",N0D,"f"),Sigma_res),"\n"))

  cat("Number of obs: ")
  cat(nrow(dat))
  cat(", groups:  ")
  nBY <- nRE <- c()
  for (i in 1:length(ranef)) { # i <- 1
    if ( (i==1) | (i>=2 & !ranef[[i]]$by %in% nBY)) {
      nBY <- c(nBY, ranef[[i]]$by)
      nRE <- c(nRE, length(unique(dat[,ranef[[i]]$by])))
    }
  }
  for (i in 1:length(nBY)) {
    if (i!=1) cat("; ")
    cat(nBY[i])
    cat(", ")
    cat(nRE[i])
  }
  cat("\n")
  #Number of obs: 60, groups:  Subj, 15; Item, 4

  cat("Fixed effects:\n")
  print(bs)

}


#==============================================================================
simLMM_verbose1 <- function(dat, Fixef_, Ranef, form_length=3) {

  Sigma_res <- Ranef[[length(Ranef)]]
  ranef <- Ranef[1:(length(Ranef)-1)]

  dat <- data.frame(dat)

  # fixed effects: obtain model matrix + write predictor names with coefficients
  if (length(Fixef_$slopes)==1 & Fixef_$slopes[1]=="1") {
    mm_fixef <- model.matrix(as.formula("~ 1"), data=dat)
  } else {
    mm_fixef <- model.matrix(as.formula(paste0("~ 1 + ", paste0(Fixef_$slopes,collapse=" + "))), data=dat)
  }
  fixNames <- colnames(mm_fixef)
  bs  <- as.character(Fixef_$b, length=form_length)
  form_fixef <- paste0("~ ")
  for (i in 1:length(bs)) {
    if (i!=1) form_fixef <- paste0(form_fixef, " + ")
    form_fixef <- paste0(form_fixef, bs[i], " * ",fixNames[i])
  }
  #lp_fixef <- mm_fixef %*% t(t(Fixef$b))


  # random effects: obtain model matrix + compute linear predictor
  form_ranef <- list()
  for (i in 1:length(ranef)) { # i <- 1
    if (length(ranef[[i]]$slopes)==1 & ranef[[i]]$slopes[1]=="1") {
      mm_ranef <- model.matrix(as.formula("~ 1"), data=dat)
    } else {
      mm_ranef <- model.matrix(as.formula(paste0("~ ",paste0(ranef[[i]]$slopes,collapse=" + "))), data=dat )
    }
    ranNames <- colnames(mm_ranef)
    SDs       <- as.character(ranef[[i]]$sd, length=form_length)
    form_ranef[[i]] <- paste0("( ")
    for (ii in 1:length(SDs)) { # ii <- 1
      if (ii!=1) form_ranef[[i]] <- paste0(form_ranef[[i]], " + ")
      form_ranef[[i]] <- paste0(form_ranef[[i]], ranNames[ii], "\\", SDs[ii])
    }
    form_ranef[[i]] <- paste0(form_ranef[[i]], " | ", ranef[[i]]$by, " )")
  }
  form_Ranef <- paste0(form_ranef, collapse=" + ")

  form_res <- paste0("1\\",as.character(Sigma_res, length=form_length))

  form_out <- paste0(form_fixef, " + ", form_Ranef, " + ", form_res, "\n")
  cat(form_out)
  for (i in 1:length(ranef)) { # i <- 1
    cat(paste0("Correlation matrix ",i,"\n"))
    print(ranef[[i]]$cor)
  }

}


#==============================================================================
simLMM_empF <- function(dat, Fixef, Ranef, family) {
#------------------------------------------------------------------------
# simulate data from a LMM using empirical=FALSE
#------------------------------------------------------------------------

Sigma_res <- Ranef[[length(Ranef)]]
if (family=="gaussian") ranef <- Ranef[1:(length(Ranef)-1)]
if (family!="gaussian") ranef <- Ranef

dat <- data.frame(dat)

# fixed effects: obtain model matrix + compute linear predictor
if (length(Fixef$slopes)==1 & Fixef$slopes[1]=="1") {
	mm_fixef <- model.matrix(as.formula("~ 1"), data=dat)
} else {
	mm_fixef <- model.matrix(as.formula(paste0("~ ", paste0(Fixef$slopes,collapse=" + "))), data=dat)
}
lp_fixef <- mm_fixef %*% t(t(Fixef$b))

# random effects: obtain model matrix + compute linear predictor
lp_ranef <- matrix(NA,nrow=nrow(dat),ncol=length(ranef))
for (i in 1:length(ranef)) { # i <- 1
	if (length(ranef[[i]]$slopes)==1 & ranef[[i]]$slopes[1]=="1") {
		mm_ranef <- model.matrix(as.formula("~ 1"), data=dat)
	} else {
		mm_ranef <- model.matrix(as.formula(paste0("~ ",paste0(ranef[[i]]$slopes,collapse=" + "))), data=dat )
	}
	# compute covariance matrix
	if (length(ranef[[i]]$sd)==1) {
		ranef[[i]]$cov <- ranef[[i]]$sd %*% ranef[[i]]$cor %*% ranef[[i]]$sd
	} else {
		ranef[[i]]$cov <- diag(ranef[[i]]$sd) %*% ranef[[i]]$cor %*% diag(ranef[[i]]$sd)
	}
	instances <- as.character(unique(dat[,ranef[[i]]$by]))
	w <- MASS::mvrnorm(length(instances), rep(0,length(ranef[[i]]$sd)), ranef[[i]]$cov, empirical=FALSE)
	for (j in 1:length(instances)) { # j <- 1
		idx <- dat[,ranef[[i]]$by]==instances[j]
		lp_ranef[idx,i] <- mm_ranef[idx,] %*% t(t(w[j,]))
	}
}
# combine linear predictions
lp <- rowSums(cbind(lp_fixef,lp_ranef))

# LINK function
#link <- "logit"
#ff <- make.link(link)
#lp <- -2:2
#ff$linkinv(as.numeric(lp))

# simulate response
if (family=="gaussian")
  resp <- lp + MASS::mvrnorm(nrow(dat), 0, Sigma_res^2, empirical=FALSE)
if (family=="lognormal")
  resp <- exp(lp + MASS::mvrnorm(nrow(dat), 0, Sigma_res^2, empirical=FALSE))
if (family=="binomial")
  resp <- as.numeric(1/(1+exp(-lp)) > runif(length(lp)))
if (family=="lp")
  resp <- as.numeric(lp)

as.numeric(resp)
}
#------------------------------------------------------



#==============================================================================
simLMM_empT <- function(dat, Fixef, Ranef, family) {
#------------------------------------------------------------------------
# simulate data from a LMM using empirical=TRUE
#------------------------------------------------------------------------

Sigma_res <- Ranef[[length(Ranef)]]
if (family=="gaussian") ranef <- Ranef[1:(length(Ranef)-1)]
if (family!="gaussian") ranef <- Ranef

dat <- data.frame(dat)

# fixed effects: obtain model matrix + compute linear predictor
mm_fixef <- model.matrix(as.formula(paste0("~ ", paste0(Fixef$slopes,collapse=" + "))), data=dat)
lp_fixef <- mm_fixef %*% t(t(Fixef$b))

# random effects: obtain model matrix + compute linear predictor
lp_ranef <- matrix(NA,nrow=nrow(dat),ncol=length(ranef))
for (i in 1:length(ranef)) { # i <- 1
	mm_ranef <- model.matrix(as.formula(paste0("~ ",paste0(ranef[[i]]$slopes,collapse=" + "))), data=dat )
	idx3 <- 1:length(Fixef$slopes)
	if (length(grep(":",Fixef$slopes))>0) idx3 <- idx3[-grep(":",Fixef$slopes)]
	if (any(Fixef$slopes=="0")) idx3 <- idx3[-which(Fixef$slopes=="0")]
	if (any(Fixef$slopes=="1")) idx3 <- idx3[-which(Fixef$slopes=="1")]
	#------------------------------------------------------------------
	if (length(idx3)>0) {
	dTmp <- data.frame(matrix(NA,nrow=nrow(dat), ncol=length(idx3)))
	for (m in 1:length(idx3)) # m <- 1
		dTmp[,m] <- as.character(factor(dat[,Fixef$slopes[idx3[m]]]))
	cell <- c()
	for (l in 1:nrow(dat)) # l <- 1
		cell[l] <- paste0(dTmp[l,], collapse=".")
	} else {
    cell <- rep(1,nrow(dat))
	}
	if (min(table(cell))<(length(ranef[[i]]$sd)+1))
	  stop('Empirical=TRUE is not possible if cell size is too small.')
	#for (l in 1:nrow(dat)) { # l <- 1
	#	cell[l] <- paste0(dat[l,Fixef$slopes[idx3]], collapse=".")
	#cell  <- paste0(dat[,Fixef$slopes[idx3]],sep=".")
	#------------------------------------------------------------------
	cellu <- unique(cell)
	table(cell)
	besselFact <- (nrow(dat)-1)/(nrow(dat)-length(unique(cellu)))
	for (k in 1:length(cellu)) { # k <- 1
		idx2 <- cell==cellu[k]
		instances <- as.character(unique(dat[idx2,ranef[[i]]$by]))
		# compute covariance matrix
		if (length(ranef[[i]]$sd)==1) {
			ranef[[i]]$cov <- ranef[[i]]$sd*sqrt(besselFact) %*% ranef[[i]]$cor %*% ranef[[i]]$sd*sqrt(besselFact)
		} else {
			ranef[[i]]$cov <- diag(ranef[[i]]$sd*sqrt(besselFact)) %*% ranef[[i]]$cor %*% diag(ranef[[i]]$sd*sqrt(besselFact))
		}
		w <- MASS::mvrnorm(length(instances), rep(0,length(ranef[[i]]$sd)), ranef[[i]]$cov, empirical=TRUE)
		for (j in 1:length(instances)) { # j <- 1
			idx <- dat[,ranef[[i]]$by]==instances[j]
			lp_ranef[idx,i] <- mm_ranef[idx,] %*% t(t(w[j,]))
		}
	}
}
# combine linear predictions
lp <- rowSums(cbind(lp_fixef,lp_ranef))

# simulate response
if (family=="gaussian" | family="lognormal") {
resp <- lp
idx3 <- 1:length(Fixef$slopes)
if (length(grep(":",Fixef$slopes))>0) idx3 <- idx3[-grep(":",Fixef$slopes)]
if (any(Fixef$slopes=="0")) idx3 <- idx3[-which(Fixef$slopes=="0")]
if (any(Fixef$slopes=="1")) idx3 <- idx3[-which(Fixef$slopes=="1")]

if (length(idx3)>0) {
dTmp <- data.frame(matrix(NA,nrow=nrow(dat), ncol=length(idx3)))
for (m in 1:length(idx3)) # m <- 1
	dTmp[,m] <- as.character(factor(dat[,Fixef$slopes[idx3[m]]]))
cell <- c()
for (l in 1:nrow(dat)) # l <- 1
	cell[l] <- paste0(dTmp[l,], collapse=".")
} else {
  cell <- rep(1,nrow(dat))
}
#if (min(table(cell))<3)
#  stop('Empirical=TRUE is not possible if cell size is too small.')
#cell <- c()
#for (l in 1:nrow(dat)) # l <- 1
#	cell[l] <- paste0(dat[l,Fixef$slopes[idx3]], collapse=".")
#cell <- paste0(dat[,Fixef$slopes],sep=".")
cellu <- unique(cell)
besselFact <- (nrow(dat)-1)/(nrow(dat)-length(unique(cellu)))
for (i in 1:length(cellu)) { # i <- 1
	idx <- cell==cellu[i]
	resp[idx] <- resp[idx] + MASS::mvrnorm(sum(idx), 0, Sigma_res^2*besselFact, empirical=TRUE)
}
}

if (family="lognormal")
  resp <- exp(resp)
if (family=="binomial")
  resp <- as.numeric(1/(1+exp(-lp)) > runif(length(lp)))
if (family=="lp")
  resp <- as.numeric(lp)


as.numeric(resp)
}
#------------------------------------------------------


#==============================================================================
extractLMM <- function(LMM) {
#------------------------------------------------------------------------
# Extract parameters (fixed and random effects) and model terms
# from a fitted linear mixed-effects model for use in model simulations
# using simLMM()
#------------------------------------------------------------------------

# extract formula from LMM
form <- as.character(LMM@call)[2]
data <- data.frame(LMM@frame)

# extract terms from formula
Terms <- extractTerms(form, data=data)

# collect fixed effects terms (+ fixed-effects estimates)
Fixef <- list(slopes=Terms[[1]], b=fixef(LMM))

# obtain terms in random effects
Ranef <- list()
for (i in 1:length(Terms[[2]])) # i <- 1
	Ranef[[i]] <- list(slopes=Terms[[2]][[i]])

# Store by, RE sd, and RE correlations
VC <- VarCorr(LMM)
nvc <- length(VC)
#idx <- list(); for (j in 1:length(Terms[[3]])) # j <- 1
#	idx[[j]] <- grep(Terms[[3]],names(VC)[j])
for (i in 1:nvc) { # i <- 1
	Ranef[[i]]$by  <- Terms[[3]][i] # names(VC)[i]
	Ranef[[i]]$sd  <- attr(VC[[i]],"stddev")
	Ranef[[i]]$cor <- attr(VC[[i]],"correlation")
}
Ranef[[nvc+1]] <- attr(VarCorr(LMM),"sc")

# return
return(Terms <- list(Fixef, Ranef, data))
}


#==============================================================================
extractTerms <- function(form, data) {
#------------------------------------------------------------------------
# Extract model terms from a linear mixed-effects formula
#------------------------------------------------------------------------

# extract terms from formula
termlab <- attr(terms(as.formula(form)),"term.labels")
interc <- attr(terms(as.formula(form)),"intercept")
# determine random effects
idx <- c(); for (i in 1:length(termlab)) # i <- 1
	idx[i] <- any(strsplit(termlab[i],"")[[1]] == "|")
# collect fixed effects terms
Fixterm <- termlab[!idx]
if (interc==1) Fixterm <- c("1", Fixterm)
if (interc==0) Fixterm <- c("0", Fixterm)
# obtain terms in random effects
ranterms <- which(idx)
#------ expand formula ------
Xterms <- c(Fixterm)
for (i in 1:length(ranterms)) { # i <- 1
	formRE <- termlab[ranterms[i]]
	tmp <- which(strsplit(formRE,"")[[1]] == "|")
	if (length(tmp)==2) {
		formRE1 <- substr(formRE, 1, tmp[1]-1)
		xterms <- c()
		# obtain information whether intercept is present
		mm <- model.matrix(as.formula(paste0("~",formRE1)), data=data)
		intercept <- apply(mm,2,function(x) all(x==1))
		if (any(intercept))
			xterms <- paste0("1", substr(formRE, tmp[length(tmp)], nchar(formRE)))
		# collect random slope terms
		REslopes <- attr(terms(as.formula(paste0("~",formRE1))),"term.labels")
		xterms <- c(xterms, paste0("0 + ",REslopes,"|", substr(formRE, tmp[length(tmp)]+1, nchar(formRE))))
		Xterms <- c(Xterms, xterms)
	} else {
		Xterms <- c(Xterms, formRE)
	}
}
termlab <- Xterms
# determine random effects
idx <- c(); for (i in 1:length(termlab)) # i <- 1
	idx[i] <- any(strsplit(termlab[i],"")[[1]] == "|")
ranterms <- which(idx)
#----------------------------
Ranterm <- list()
BY <- c()
for (i in 1:length(ranterms)) { # i <- 1
	formRE <- termlab[ranterms[i]]
	# extract random slopes
	formRE1 <- substr(formRE,1,which(strsplit(formRE,"")[[1]] == "|")[1]-1)
	# obtain information whether intercept is present
	mm <- model.matrix(as.formula(paste0("~",formRE1)), data=data)
	intercept <- apply(mm,2,function(x) all(x==1))
	if (any(intercept)) {intercept <- "1"} else {intercept <- "0"}
	# collect random slope terms
	Ranterm[[i]] <- c(intercept, attr(terms(as.formula(paste0("~",formRE1))),"term.labels"))
	if (length(Ranterm[[i]])==0) Ranterm[[i]] <- "1"
	# extract BY argument
	tmp <- which(strsplit(formRE,"")[[1]] == "|")
	formRE2 <- substr(formRE,tmp[length(tmp)]+1,nchar(formRE))
	BY[i]   <- attr(terms(as.formula(paste0("~",formRE2))),"term.labels")
}

(Terms <- list(Fixterm, Ranterm, BY))
}


