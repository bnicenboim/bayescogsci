library(MASS)

genfake <- function(dat=data,
                    nsubj=37,nitem=15,
                    beta0=5.88,beta1=0.06,
                    sigma_u0=.33,
                    sigma_u1=.11,
                    sigma_w0=.04,
                    sigma_w1=.04,
                    rho_u=0.51,
                    rho_w=.29,
                    sigma=.2,
                    N=dim(dat)[1]) {
  Sigma_u <- matrix(c(sigma_u0^2, sigma_u0*sigma_u1*rho_u,
                      sigma_u0*sigma_u1*rho_u,sigma_u1^2),
                    ncol=2)
  ## generate by subject intercepts and slopes:
  u <- mvrnorm(nsubj,c(0,0),Sigma_u)
  Sigma_w <- matrix(c(sigma_w0^2, sigma_w0*sigma_w1*rho_w,
                      sigma_w0*sigma_w1*rho_w,sigma_w1^2),
                    ncol=2)
  ## generate by item intercepts and slopes:
  w <- mvrnorm(nitem,c(0,0),Sigma_w)
  ## generate normally distributed data:
  rtfake <- NA
  for (j in 1:N){ # j <- 1
    rtfake[j] <- rlnorm(1, beta0 + 
                          u[dat$subj[j],1] +
                          w[dat$item[j],1] +
                          (beta1 +
                             u[dat$subj[j],2] +
                             w[dat$item[j],2]
                          ) * dat$x[j],
                        sigma)
  }
  rtfake
}
