#' Implement the RPSFT + Cox method
#'
#' This function implements the proposed RPSFT + Cox method and returns an 
#' estimate of the treatment effect on the overall survival. It can be used if
#' investigators have sufficient knowledge of the PDT effects and the survival 
#' time during the PDT period.
#' 
#'
#' @param t.pfs the progression-free survival (PFS) time
#' @param t.co the observed time from the start to the end of the crossover 
#' period
#' @param t.os the overall survival (OS) time
#' @param delta.os the event indicator of the OS
#' @param cen.time the censoring time. Under administrative censoring, the 
#' censoring time is available for each subject
#' @param a the initial binary treatment indicator encoded as 0 (control) or 
#' 1 (treatment)
#' @param mat.pdt the matrix of all prognosis factors that may affect the PDTs, 
#' dimension: n*p, n: the sample size, p: dimension of the prognosis factors
#' @param delta.co the vector of crossover indicators
#' @param delta.pdt the vector of the PDT indicators
#' @param include.pdt a logical variable to indicate the inclusion of the PDT
#' effects in the OS analysis. Default: `include.pdt = TRUE`.
#' @param const the multiplicative sensitivity parameter that reflects a change
#' in the treatment benefit for subjects in the control group with treatment 
#' crossover.
#' @param tau.lower the lower bound of the treatment effect parameter for the 
#' grid search. Default: `tau.lower = -1`.
#' @param tau.upper the upper bound of the treatment effect parameter for the 
#' grid search. Default: `tau.upper = 1`.
#' @param grid.length the number of values in the grid search. Default: 
#' `grid.length = 50`.
#' @return tau.est an estimate of the treatment effect
#' @import survival
#' @export
#' 
#' @examples
#' x <- cbind(dat$x1, dat$x2)
#' a <- dat$a
#' t.pfs <- dat$t.pfs
#' t.co <- dat$t.co
#' t.os <- dat$t.os
#' cen.time <- dat$c
#' delta.pfs <- dat$delta.pfs
#' delta.dp <- dat$delta.dp
#' delta.pdt <- dat$delta.pdt
#' delta.co <- dat$delta.co
#' delta.os <- dat$delta.os
#' mat.pdt <- x 
#' res.rpsftcox <- rpsft.cox(t.pfs, t.co, t.os, delta.os, 
#'                           cen.time, a, mat.pdt, delta.co, delta.pdt,
#'                           include.pdt = TRUE, const = 1,
#'                           tau.lower = -3, tau.upper = 1,
#'                           grid.length = 100)
#' res.rpsftcox
rpsft.cox <- function(t.pfs, t.co, t.os, delta.os, 
                      cen.time, a, mat.pdt, delta.co, delta.pdt,
                      include.pdt = TRUE, const = 1,
                      tau.lower = -1, tau.upper = 1,
                      grid.length = 50, var.boot = FALSE, B = NULL){
  
  # dimension of prognosis factors
  p <- ncol(mat.pdt)
  # length of the data
  n <- length(t.pfs)
  
  if(include.pdt){
    ## Step 1: estimate the PDT effect
    index.pdt <- which(t.co < t.os)
    
    t.pdp <- (t.os - t.co)[index.pdt]
    dat.pdt <- data.frame(t.pdt = t.pdp, delta.os = delta.os[index.pdt],
                          a = a[index.pdt], delta.co = delta.co[index.pdt],
                          aco = a[index.pdt] + delta.co[index.pdt],
                          delta.pdt = delta.pdt[index.pdt])
    mat0.pdt <- cbind(1, mat.pdt[index.pdt,])
    mat0.pdt.int <- apply(mat0.pdt, 2, function(x) dat.pdt$delta.pdt*x)
    # Step 1: estimate the PDT effects
    # fit a stratified cox model
    surv.obj <- with(dat.pdt, survival::Surv(t.pdt, delta.os))
    fit.cox <- survival::coxph(surv.obj ~ mat0.pdt.int +
                                 survival::strata(aco), 
                               data = dat.pdt)
    linear.pred <- predict(fit.cox, newdata = dat.pdt)
    pdtmodel.long <- rep(0, n)
    pdtmodel.long[index.pdt] <- linear.pred
  }
  
  else{
    pdtmodel.long <- rep(0, n)
  }
  
  
  ## Step 2: estimate the treatment effect
  res.final <- rpsft.linear(a = a, t.pfs = t.pfs, t.co = t.co, t.os = t.os,
                            delta.os = delta.os, delta.co = delta.co,
                            cen.time = cen.time,
                            pdt.effect = pdtmodel.long,
                            const = const,
                            tau.lower = tau.lower, tau.upper = tau.upper,
                            grid.length = grid.length)
  tau.final <- res.final$tau
  
  ## Step 3: 
  res.list <- list(tau.est = tau.final)
  
  return(list(
    tau.est = tau.final
  ))
}

#' Implement the RPSFT + IPCW method
#'
#' This function implements the proposed RPSFT + IPCW method and returns an 
#' estimate of the treatment effect on the overall survival. It can be used if
#' investigators have sufficient knowledge of the confounding factors that may
#' affect the probability of receiving the PDTs.
#'
#' @param t.pfs the progression-free survival (PFS) time
#' @param t.co the observed time from the start to the end of the crossover 
#' period
#' @param t.os the overall survival (OS) time
#' @param delta.os the event indicator of the OS
#' @param cen.time the censoring time. Under administrative censoring, the 
#' censoring time is available for each subject
#' @param a the initial binary treatment indicator encoded as 0 (control) or 
#' 1 (treatment)
#' @param mat.pdt the matrix of all prognosis factors that may affect the 
#' patient to receive the PDTs, with the dimension as n*p, where n: the sample 
#' size, p: dimension of the prognosis factors. Inclusion of the survival times
#' before the PDT period is recommended.
#' @param delta.co the event indicator of treatment crossover
#' @param delta.pdt the event indicator of the PDTs
#' @param include.pdt a logical variable to indicate the inclusion of the PDT
#' effects in the OS analysis. Default: `include.pdt = TRUE`.
#' @param const the multiplicative sensitivity parameter that reflects a change
#' in the treatment benefit for subjects in the control group with treatment 
#' crossover.
#' @param tau.lower the lower bound of the treatment effect parameter for the 
#' grid search. Default: `tau.lower = -1`.
#' @param tau.upper the upper bound of the treatment effect parameter for the 
#' grid search. Default: `tau.upper = 1`.
#' @param grid.length the number of values in the grid search. Default: 
#' `grid.length = 50`.
#' @return tau.est an estimate of the treatment effect
#' @export
#' @examples
#' x <- cbind(dat$x1, dat$x2)
#' a <- dat$a
#' t.pfs <- dat$t.pfs
#' t.co <- dat$t.co
#' t.os <- dat$t.os
#' cen.time <- dat$c
#' delta.pfs <- dat$delta.pfs
#' delta.dp <- dat$delta.dp
#' delta.pdt <- dat$delta.pdt
#' delta.co <- dat$delta.co
#' delta.os <- dat$delta.os
#' mat.pdt <- cbind(x, t.co)
#' res.rpsftipcw <- rpsft.ipcw(t.pfs, t.co, t.os, delta.os, 
#'                           cen.time, a, mat.pdt, delta.co, delta.pdt,
#'                           include.pdt = TRUE, const = 1,
#'                           tau.lower = -3, tau.upper = 1,
#'                           grid.length = 100)
#' res.rpsftipcw
rpsft.ipcw <- function(t.pfs, t.co, t.os, delta.os, 
                       cen.time, a, mat.pdt, delta.co, delta.pdt,
                       include.pdt = TRUE, const = 1,
                       tau.lower = -1, tau.upper = 1,
                       grid.length = 50){
  # dimension of prognosis factors
  p <- ncol(mat.pdt)
  # length of the data
  n <- length(t.pfs)
  
  # transform mat.pdt to a data.frame
  df.mat.pdt <- data.frame(mat.pdt)
  colnames(df.mat.pdt) <- paste0("x", 1:p) 
  
  ## Step 1: compute IPCW to account for the PDT effect
  if(include.pdt){
    ## for those who experience the PDTs
    index.pdt <- which(t.co < t.os)
    dat.pdt <- data.frame(a = a[index.pdt], delta.co = delta.co[index.pdt],
                          delta.pdt = delta.pdt[index.pdt])
    dat0.pdt <- cbind(dat.pdt, df.mat.pdt[index.pdt,])
    mat0.pdt <- df.mat.pdt[index.pdt,]
    
    ## for those who would have experienced the PDTs (those who progressed)
    index.pdt0 <- which(t.pfs < t.os)
    dat.pdt0 <- data.frame(a = a[index.pdt0], delta.co = delta.co[index.pdt0],
                           delta.pdt = delta.pdt[index.pdt0])
    dat0.pdt0 <- cbind(dat.pdt0, df.mat.pdt[index.pdt0,])
    
    mat0.pdt0 <- df.mat.pdt[index.pdt0,]
    
    formula.logit <- reformulate(paste0("x", 1:p), "delta.pdt")
    # (1) A = 1
    fit.pdt.a1 <- glm(formula.logit,
                      family = binomial,
                      data = dat0.pdt[dat0.pdt$a == 1,])
    prob.a1 <- predict(fit.pdt.a1, newdata = dat0.pdt0[dat0.pdt0$a == 1,],
                       type = "response")
    # (2) A = 0 & delta.co = 0
    fit.pdt.a0noco <- glm(formula.logit,
                          family = binomial,
                          data = dat0.pdt[dat0.pdt$a == 0 & 
                                            dat0.pdt$delta.co == 0,])
    prob.a0noco <- predict(fit.pdt.a0noco, 
                           newdata = dat0.pdt0[dat0.pdt0$a == 0 & 
                                                 dat0.pdt0$delta.co == 0,],
                           type = "response")
    # (3) A = 0 & delta.co = 1
    fit.pdt.a0co <- glm(formula.logit,
                        family = binomial,
                        data = dat0.pdt[dat0.pdt$a == 0 & 
                                          dat0.pdt$delta.co == 1,])
    prob.a0co <- predict(fit.pdt.a0co, 
                         newdata = dat0.pdt0[dat0.pdt0$a == 0 & 
                                               dat0.pdt0$delta.co == 1,],
                         type = "response")
    
    # (1) unstablized weights
    ipcw.comb <- rep(1, n)
    ipcw.comb[a == 1 & t.pfs < t.os] <- 1/(1 - prob.a1)
    ipcw.comb[a == 0 & delta.co == 0 & t.pfs < t.os] <- 1/(1 - prob.a0noco)
    ipcw.comb[a == 0 & delta.co == 1 & t.pfs < t.os] <- 1/(1 - prob.a0co)
    
    # (2) stablized weights (need further check)
    # 
    # ipcw.comb <- rep(1, n)
    # ipcw.comb[a == 1 & t.co < t.os] <- 1/(1 - prob.a1)
    # ipcw.comb[a == 0 & delta.co == 0 & t.co < t.os] <- 1/(1 - prob.a0noco)
    # ipcw.comb[a == 0 & delta.co == 1 & t.co < t.os] <- 1/(1 - prob.a0co)
  }
  
  else{
    ipcw.comb <- rep(1,n)
  }
  
  ## Step 2: estimate the treatment effect
  res.final <- rpsft.weighted(a = a, t.pfs = t.pfs, t.co = t.co, t.os = t.os,
                              delta.os = delta.os, delta.co = delta.co,
                              delta.pdt = delta.pdt, cen.time = cen.time,
                              const = 1,
                              weights = ipcw.comb,
                              tau.lower = tau.lower, tau.upper = tau.upper,
                              grid.length = grid.length)
  tau.final <- res.final$tau
  
  ## Step 3: estimate the variance (bootstrap)
  
  
  
  return(list(
    tau.est = tau.final
  ))
}

#' Fit the RPSFT model on the original data with the estimated PDT effects
#'
#' This function fits the RPSFT model and returns the treatment effect estimate 
#' and the corresponding p-value which result in a zero test statistic of the 
#' log-rank test. It is an intermediate function that corresponds to the second 
#' step of the proposed RPSFT + Cox method to account for the treatment
#' crossover. It can be used if investigators have sufficient knowledge of the 
#' PDT effects and the survival time during the PDT period.
#'
#' @param t.pfs the progression-free survival (PFS) time
#' @param t.co the observed time from the start to the end of the crossover 
#' period
#' @param t.os the overall survival (OS) time
#' @param a the initial binary treatment indicator encoded as 0 (control) or 
#' 1 (treatment)
#' @param delta.co the event indicator of treatment crossover
#' @param delta.os the event indicator of the OS
#' @param cen.time the censoring time. Under administrative censoring, the 
#' censoring time is available for each subject
#' @param pdt.effect a vector of the estimated PDT effects
#' @param const the multiplicative sensitivity parameter that reflects a change
#' in the treatment benefit for subjects in the control group with treatment 
#' crossover. Default: `const = 1`.
#' @param weights the use-defined weights used in the log-rank test. Default: 
#' each individual has an equal weight, i.e., `weights = rep(1, length(a))`
#' @param tau.lower the lower bound of the treatment effect parameter for the 
#' grid search. Default: `tau.lower = -1`.
#' @param tau.upper the upper bound of the treatment effect parameter for the 
#' grid search. Default: `tau.upper = 1`.
#' @param grid.length the number of values in the grid search. Default: 
#' `grid.length = 50`.
#' @return a list of test results
#' \itemize{
#' \item test.stat: the test statistic
#' \item p.value: the p-value
#' \item tau: the estimate that results in a zero test statistic
#' \item test.stat.seq: a sequence of test statistics in the grid search
#' \item p.value.seq: a sequence of p-values in the grid search
#' }
#' @export
rpsft.linear <- function(t.pfs, t.co, t.os, a,
                         delta.os, delta.co, cen.time,
                         pdt.effect, const = 1,
                         weights = rep(1, length(a)),
                         tau.lower = -1, tau.upper = 1, grid.length = 50){
  t.pdt.a0co.w <- (t.os - t.co)*exp(pdt.effect)
  t.pdt.noco.w <- (t.os - t.pfs)*exp(pdt.effect)
  t.pdp.a0co.w <- t.co - t.pfs + t.pdt.a0co.w
  
  grid.tau <- seq(tau.lower, tau.upper, length = grid.length)
  test.stat <- c(0)
  p.value <- c(0)
  for(i in 1:grid.length){
    tau.temp <- grid.tau[i]
    u.a1 <- (t.pfs + t.pdt.noco.w)*exp(tau.temp)
    u.a0co <- t.pfs + t.pdp.a0co.w*exp(const*tau.temp)
    u.a0nco <- t.pfs + t.pdt.noco.w
    
    c.adj1 <- (t.pfs+ (cen.time - t.pfs)*exp(pdt.effect))*exp(tau.temp)
    c.adj2 <- cen.time*exp(tau.temp)
    c.adj3 <- cen.time
    c.adj4 <- t.pfs+ (cen.time - t.pfs)*exp(pdt.effect)
    c.adj5 <- t.pfs + (cen.time - t.pfs)*exp(const*tau.temp)
    c.adj6 <- t.pfs + (t.co - t.pfs + (cen.time - t.co)*exp(pdt.effect))*
      exp(const*tau.temp)
    
    c.adj <- apply(cbind(c.adj1, c.adj2, c.adj3, c.adj4, c.adj5, c.adj6), 1, 
                   min)
    
    u.adj <- (1-a)*(delta.co*u.a0co + (1 - delta.co)*u.a0nco) + a*u.a1
    uc.adj <- cbind(u.adj, c.adj)
    v1.adj <- apply(uc.adj, 1, function(x) min(c(x[1],x[2])))
    delta.os.adj <- delta.os*apply(uc.adj, 1, function(x) as.numeric(x[1]<x[2]))
    t.os.adj <- (1 - delta.os.adj)*c.adj + delta.os.adj*v1.adj
    order.ind <- sort(t.os.adj, index.return = TRUE)$ix
    unordered.weight <- weights
    order.weight <- unordered.weight[order.ind]
    
    data.adj <- data.frame(t.os = t.os.adj, 
                           delta.os = delta.os.adj, 
                           a = a)
    test.res <- weighted.logrank(Surv(t.os, delta.os) ~ a, 
                                 data = data.adj,
                                 weights = order.weight)
    test.stat[i] <- test.res$test.stat
    p.value[i] <- test.res$p.value
  }
  index.min <- which.min(test.stat)
  test.stat.min <- test.stat[index.min]
  p.value.opt <- p.value[index.min]
  tau.opt <- grid.tau[index.min]
  
  return(list(test.stat = test.stat.min,
              p.value = p.value.opt,
              tau = tau.opt,
              test.stat.seq = test.stat,
              p.value.seq = p.value))
}

#' Fit the RPSFT model on the original data with the estimated PDT effects
#'
#' This function fits the RPSFT model and returns the treatment effect estimate 
#' and the corresponding p-value which result in a zero test statistic of the 
#' log-rank test. It is an intermediate function that corresponds to the second 
#' step of the proposed RPSFT + Cox method to account for the treatment
#' crossover. It can be used if investigators have sufficient knowledge of the 
#' PDT effects and the survival time during the PDT period.
#'
#' @param t.pfs the progression-free survival (PFS) time
#' @param t.co the observed time from the start to the end of the crossover 
#' period
#' @param t.os the overall survival (OS) time
#' @param a the initial binary treatment indicator encoded as 0 (control) or 
#' 1 (treatment)
#' @param delta.co the event indicator of treatment crossover
#' @param delta.pdt the event indicator of the PDTs
#' @param delta.os the event indicator of the OS
#' @param cen.time the censoring time. Under administrative censoring, the 
#' censoring time is available for each subject
#' @param const the multiplicative sensitivity parameter that reflects a change
#' in the treatment benefit for subjects in the control group with treatment 
#' crossover. Default: `const = 1`.
#' @param weights the use-defined weights used in the log-rank test. Default: 
#' each individual has an equal weight, i.e., `weights = rep(1, length(a))`
#' @param tau.lower the lower bound of the treatment effect parameter for the 
#' grid search. Default: `tau.lower = -1`.
#' @param tau.upper the upper bound of the treatment effect parameter for the 
#' grid search. Default: `tau.upper = 1`.
#' @param grid.length the number of values in the grid search. Default: 
#' `grid.length = 50`.
#' @return a list of test results
#' \itemize{
#' \item test.stat: the test statistic
#' \item p.value: the p-value
#' \item tau: the estimate that results in a zero test statistic
#' \item test.stat.seq: a sequence of test statistics in the grid search
#' \item p.value.seq: a sequence of p-values in the grid search
#' }
#' @export
rpsft.weighted <- function(t.pfs, t.co, t.os, a,
                           delta.os, delta.co, delta.pdt, cen.time,
                           const = 1,
                           weights = rep(1, length(a)),
                           tau.lower = -1, tau.upper = 1, grid.length = 50){
  t.os.adj0 <- delta.pdt*t.co + (1 - delta.pdt)*t.os
  t.pdp <- t.os.adj0 - t.pfs
  delta.os.adj0 <- (1 - delta.pdt)*delta.os
  
  grid.tau <- seq(tau.lower, tau.upper, length = grid.length)
  test.stat <- c(0)
  p.value <- c(0)
  for(i in 1:grid.length){
    tau.temp <- grid.tau[i]
    u.a1 <- (t.pfs + t.pdp)*exp(tau.temp)
    u.a0co <- t.pfs + t.pdp*exp(const*tau.temp)
    u.a0nco <- t.pfs + t.pdp
    
    c.adj1 <- cen.time*exp(tau.temp)
    c.adj2 <- cen.time
    c.adj3 <- t.pfs + (cen.time - t.pfs)*exp(const*tau.temp)
    c.adj <- apply(cbind(c.adj1, c.adj2, c.adj3), 1, min)
    
    u.adj <- (1-a)*(delta.co*u.a0co + (1 - delta.co)*u.a0nco) + a*u.a1
    uc.adj <- cbind(u.adj, c.adj)
    v1.adj <- apply(uc.adj, 1, function(x) min(c(x[1],x[2])))
    delta.os.adj <- delta.os.adj0*apply(uc.adj, 1, function(x) as.numeric(x[1]<x[2]))
    t.os.adj <- (1 - delta.os.adj)*c.adj + delta.os.adj*v1.adj
    order.ind <- sort(t.os.adj, index.return = TRUE)$ix
    unordered.weight <- weights
    order.weight <- unordered.weight[order.ind]
    
    data.adj <- data.frame(t.os = t.os.adj, 
                           delta.os = delta.os.adj, 
                           a = a)
    test.res <- weighted.logrank(Surv(t.os, delta.os) ~ a, 
                                 data = data.adj,
                                 weights = order.weight)
    test.stat[i] <- test.res$test.stat
    p.value[i] <- test.res$p.value
  }
  index.min <- which.min(test.stat)
  test.stat.min <- test.stat[index.min]
  p.value.opt <- p.value[index.min]
  tau.opt <- grid.tau[index.min]
  
  return(list(test.stat = test.stat.min,
              p.value = p.value.opt,
              tau = tau.opt,
              test.stat.seq = test.stat,
              p.value.seq = p.value))
}

#' Conduct the weighted log-rank test
#'
#' This function implements the log-rank test with user-defined weights and 
#' returns the test statistics and the p-value. 
#'
#' @param formula the formula needed for the log-rank test
#' @param data the data frame consists of the survival time, the event 
#' indicator, and the initial treatment assignment
#' @param weights the use-defined weights used in the log-rank test. Default: 
#' each individual has an equal weight, i.e., `weights = rep(1, length(a))`
#' @return the test statistic and the p-value of the weighted log-rank test
#' @import nphRCT
#' @export
weighted.logrank <- function(formula, data, weights = rep(1, nrow(data))){
  info.test <- nphRCT::find_at_risk(formula, data)
  
  n_risk1 <- info.test$n_risk[1:(nrow(info.test) - 1)]
  n_risk0 <- info.test$n_risk[2:nrow(info.test)]
  diff.num <- n_risk1 - n_risk0
  temp <- c(cumsum(diff.num), cumsum(diff.num)[length(cumsum(diff.num))]+1)
  weights.shorten <- weights[temp]
  
  Uw <- sum(with(info.test, weights.shorten*(n_event_0 - n_event*n_risk_0/n_risk)))
  Vw.all <- with(info.test, ifelse(n_risk_0*n_risk_1 == 0, 0, n_risk_0*n_risk_1*n_event*(n_risk - n_event)/(n_risk^2*(n_risk - 1))))
  Vw <- sum(weights.shorten^2*Vw.all)
  
  test.stat <- (Uw/sqrt(Vw))^2
  p.value <- 1 - pchisq(test.stat, 1)
  
  return(list(test.stat = test.stat,
              p.value = p.value))
}

#' @importFrom stats binomial glm pchisq predict reformulate
NULL