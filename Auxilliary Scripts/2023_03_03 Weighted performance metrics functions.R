#---------------------------------------------------
# Description: Functions for the computation of weighted performance metrics
# Author: B. Rentroia Pacheco
#---------------------------------------------------

# Truncate Follow-up time
#---------------------------------------------------
# @param outcome_FUP: numeric vector with follow up time of subjects
# @param outcome_num: numeric vector with status indicator, 0 =alive, 1 = dead. Follows the same conventions as parameter event in the Surv function from the survival package.
# @param tp: number indicating at which timepoint truncation should occur. It should be between 0 and max of outcome_FUP.
truncateFUP <-function (outcome_FUP,outcome_num,tp){
  outcome_num_u <- outcome_num
  i.greaterFUP <- which(outcome_FUP>tp)
  if(length(i.greaterFUP)>0){
    outcome_num_u[outcome_FUP>tp]<-0
  }
  outcome_FUP[outcome_FUP>tp] <- tp
  return(list("FUP" = outcome_FUP,"Status"= outcome_num_u))
}

# Discriminative ability measured with Harrel's C-index:
#---------------------------------------------------
# @param outcome_FUP: numeric vector with follow up time of subjects
# @param outcome_num: numeric vector with status indicator, 0 =alive, 1 = dead. Follows the same conventions as parameter event in the Surv function from the survival package.
# @param predictions: numeric vector with model risk predictions 
# @param tp: number indicating at which timepoint the OE ratio should be computed. It should be between 0 and max of outcome_FUP.
# @param weights: if required output is weighted C-index, weights is a numeric vector with weights of each observation. If required output is unweighted C-index, weights should be NULL
computeCindex <- function(outcome_FUP,outcome_num,predictions,tp,weights=NULL){
  #Truncate outcome at timepoint t:
  if (!is.null(tp)){
    truncated_survival <- truncateFUP(outcome_FUP,outcome_num,tp)
    outcome_FUP <- truncated_survival[["FUP"]]
    outcome_num <- truncated_survival[["Status"]]
  }
  
  # Compute unweighted and weighted C-index:
  if (is.null(weights)){
    cind=1-rcorr.cens(predictions, Surv(outcome_FUP,outcome_num))[1]
  }else{
    cind=intsurv::cIndex(outcome_FUP,event=outcome_num,predictions,weight=weights)[[1]]
  }
  names(cind)=NULL
  return(cind)
}

# Discriminative ability metrics: SE, SP, NPV,PPV, for survival data
#---------------------------------------------------
# @param outcome_FUP: numeric vector with follow up time of subjects
# @param outcome_num: numeric vector with status indicator, 0 =alive, 1 = dead. Follows the same conventions as parameter event in the Surv function from the survival package.
# @param predictions: numeric vector with model risk predictions 
# @param cutoff: prediction threshold
# @param tp: number indicating at which timepoint the OE ratio should be computed. It should be between 0 and max of outcome_FUP.
# @param weights: if required output is weighted metrics, weights is a numeric vector with weights of each observation. If required output is unweighted C-index, weights should be NULL
computeDiscMetrics <- function(outcome_FUP,outcome_num,predictions,cutoff,tp,weights=NULL){
  
  # Obtain labels for the desired cutoff:
  labels_cutoff <- ifelse(predictions>cutoff,1,0)
  
  # Weights are all equal if they are not specified:
  if (is.null(weights)){
    weights <-rep(1,length(predictions))
  }
  
  # Compute TP, TN, FP and FN:
  low_risk_surv <-summary(survfit(Surv(outcome_FUP[labels_cutoff ==0],outcome_num[labels_cutoff ==0])~1,weights=weights[labels_cutoff ==0]),time=tp)$surv
  high_risk_surv <-summary(survfit(Surv(outcome_FUP[labels_cutoff ==1],outcome_num[labels_cutoff ==1])~1,weights=weights[labels_cutoff ==1]),time=tp)$surv
  high_risk_n <- sum(ifelse(labels_cutoff ==1,weights,0))
  low_risk_n <-sum(ifelse(labels_cutoff ==0,weights,0))
  overall_surv <- summary(survfit(Surv(outcome_FUP,outcome_num)~1,weights=weights),time=tp)$surv

  # TN: 
  TN <-   low_risk_surv * low_risk_n
  # TP: 
  TP <-  (1-high_risk_surv ) *high_risk_n
  # FN:
  FN <-  (1-low_risk_surv) * low_risk_n
  # FP:
  FP <-  high_risk_surv  * high_risk_n
  # All:
  Pos <-  (1-overall_surv) * sum(weights)
  Neg <-  overall_surv * sum(weights)
  
  # Compute Sensitivity:
  SE <- TP/(Pos) # same output as survivalROC, method = "KM". Note we obtain similar results with the denominator = TP+FN, but not equal.
  
  # Compute Specificity:
  SP <- 1-FP/Neg  # same output as survivalROC, method = "KM". Note we obtain similar results with the formula TN/(TN+FP), but not exactly the same.
  
  # Compute PPV:
  PPV <-  (1-high_risk_surv) #as per the definiton. if you rewrite PPV = TP/(TP+FP), you will see that many terms cancel out
  # Compute NPV:
  NPV <- low_risk_surv #as per the definiton. if you rewrite NPV = TN/(TN+FN), you will see that many terms cancel out
  
  # Summarize all metrics
  performance_metrics <-c("SE"=SE,"SP"=SP,"PPV"=PPV,"NPV"=NPV)
  
  return(performance_metrics)
}


# Observed to events ratio at timepoint tp:
#---------------------------------------------------
# @param outcome_FUP: numeric vector with follow up time of subjects
# @param outcome_num: numeric vector with status indicator, 0 =alive, 1 = dead. Follows the same conventions as parameter event in the Surv function from the survival package.
# @param predictions: numeric vector with model risk predictions 
# @param tp: number indicating at which timepoint the OE ratio should be computed. It should be between 0 and max of outcome_FUP.
# @param weights: if required output is weighted OE ratio, weights is a numeric vector with weights of each observation. If required output is unweighted OE ratio, weights should be NULL
computeOEratio <- function(outcome_FUP,outcome_num,predictions,tp,weights=NULL){
  # unweighted and weighted version:
  if (is.null(weights)){
    obj <- summary(survfit(Surv(outcome_FUP,outcome_num) ~ 1), 
                   times = tp,extend=TRUE)
    OE_ratio <-  (1 - obj$surv) / (mean(predictions))
  }else{
    obj<- summary(survfit(Surv(outcome_FUP,outcome_num) ~ 1,weights = weights), 
                  times = tp,extend=TRUE)
    OE_ratio <-  (1 - obj$surv) / weighted.mean(predictions,weights)
  }
  return(OE_ratio)
}

# Calibration slope at timepoint tp:
#---------------------------------------------------
# @param outcome_FUP: numeric vector with follow up time of subjects
# @param outcome_num: numeric vector with status indicator, 0 =alive, 1 = dead. Follows the same conventions as parameter event in the Surv function from the survival package.
# @param predictions: numeric vector with model predictions
# @param tp: number indicating at which timepoint the data should be truncated at.  It should be between 0 and max of outcome_FUP.
# @param weights: if required output is weighted calibration slope, weights is a numeric vector with weights of each observation. If required output is unweighted calibration slope, weights should be NULL
calibrationSlope <- function(outcome_FUP,outcome_num,predictions,tp=NULL,weights=NULL){
    #Truncate outcome at timepoint t:
    if (!is.null(tp)){
      truncated_survival <- truncateFUP(outcome_FUP,outcome_num,tp)
      outcome_FUP <- truncated_survival[["FUP"]]
      outcome_num <- truncated_survival[["Status"]]
    }
    # Compute calibration slope:
    log_log_pred <- log(-log(1-predictions)) #Log log of the survival predictions
    cal_slope<- coef(coxph(Surv(outcome_FUP,outcome_num)~log_log_pred,weights=weights))
    names(cal_slope)<- NULL
  return(cal_slope)
}


# Calibration plot, with truncation of survival events/times at timepoint t
#---------------------------------------------------
# @param predictions: numeric vector with model risk predictions
# @param outcome_num: numeric vector with status indicator, 0 =alive, 1 = dead. Follows the same conventions as parameter event in the Surv function from the survival package.
# @param outcome_FUP: numeric vector with follow up time of subjects
# @param u: number indicating at which timepoint the data should be truncated at.  It should be between 0 and max of outcome_FUP.
# @param bin_nr: Approximate number of subjects per group
# @param lim: Vector of length 2, indicating the limits of the axes, e.g. c(0,1)
# @param pl: TRUE if user wants to plot the calibration plot, FALSE otherwise.
# @param weights: if required output is weighted calibration plot, weights is a numeric vector with weights of each observation. If required output is unweighted calibration plot, weights should be NULL
# @param log_plot: TRUE if axes should be in log scale, FALSE otherwise
# @param main_plot: Title of the plot
# @param surv_or_risk: String with two possible values: "surv" if survival probabilities should be shown, or "risk" if risk probabilities should be shown.
# @param new_plot: TRUE if new plot should be created, or FALSE if calibration plot will be added to a pre existent plot.
# @param dots_col: String indicating the color of the dots in the calibration plot
# @param show_segments: TRUE if histogram with distribution of risk/survival probabilities should be shown, FALSE otherwise
calib_plot =function(predictions,outcome_num,outcome_FUP,u,bin_nr,lim,pl,weights=NULL,log_plot=FALSE,main_plot="",surv_or_risk = "surv",new_plot = TRUE,dots_col = "darkred",show_segments = FALSE){
  
  # Check that plot parameters are compatible:
  
  
  if(show_segments & (!pl)){
    print("Segments cannot be shown if plot parameter pl is set to FALSE")
  }
  
  # Truncate outcome at timepoint t:
  if (!is.null(u)){
    truncated_survival <- truncateFUP(outcome_FUP,outcome_num,u)
    outcome_FUP <- truncated_survival[["FUP"]]
    outcome_num <- truncated_survival[["Status"]]
  }
  
  # Calculate calibration metrics:
  fit_intercept <- computeOEratio(outcome_FUP,outcome_num,predictions,u,weights) 
  fit_slope<-calibrationSlope(outcome_FUP,outcome_num,predictions,NULL,weights) 
  
  # Axes names:
  if(surv_or_risk =="risk"){
    ylab <- paste0("1-Kaplan-Meier survival ",u,' years') 
    xlab <- "Predicted risk probability"
  }else if (surv_or_risk =="surv"){
    ylab <- paste0("Kaplan-Meier survival ",u,' years') 
    xlab <- "Predicted survival probability"
  }
  
  # Plot:
  
  if(new_plot){
    if(log_plot){
      # Axes limits cannot contain 0
      if (lim[1]==0){
        print("Warning: Axes limits cannot contain 0 in log plots. Lower limit was replaced by 0.0001")
        lim[1] <-0.0001
      }
      plot(.5, .5, xlim=lim, ylim=lim, type="n", xlab=xlab, ylab=ylab,log ="xy",main=main_plot,cex.lab=1.5, cex.axis=1.5)
    }else{
      plot(.5, .5, xlim=lim, ylim=lim, type="n", xlab=xlab, ylab=ylab,main=main_plot,cex.lab=1.5, cex.axis=1.5)
    }
    abline(0, 1, lwd=6, col=gray(.85))
  }
  
  
  
  lt <- 1; leg <- "Ideal"; marks <- -1; lwd <- 6; col <- gray(.85)
  
  # Partition data, to show groups in the calibration plot:
  # Code was adapted from rms's function groupKM, this could not be used directly because it didn't have weighted KM estimates
  
  # Predictions should be survival predictions:
  predictions = 1-predictions
  
  # Divide observations:
  
  if(is.null(weights)){
    bin <- cut2(predictions,m=bin_nr,levels.mean=TRUE,digits=7)
  }else{
    nr_percentiles <-floor(length(predictions)/bin_nr)
    bin <- cut2(seq(0.0001,1,0.0001),m=10000/nr_percentiles,levels.mean=TRUE,digits=2) #obtain the desired percentiles
    pct_bin <-round(as.numeric(levels(bin))+as.numeric(levels(bin))[1],2)
    upper_bond_qnt <- as.numeric(wtd.quantile(predictions,w=weights,probs=c(pct_bin)))
    bin<-cut2(predictions,c(0,upper_bond_qnt),levels.mean=TRUE)
  } 
  # Setup:
  Srv = Surv(outcome_FUP,outcome_num)
  means <- as.numeric(levels(bin))
  q <- unclass(bin)  # to assign observations to groups
  g <- length(levels(q))
  
  km      <- double(g)
  pred    <- km
  std.err <- km
  events  <- integer(g)
  numobs  <- events
  ses_mean <- double(g)
  
  #Compute KM estimates for each group:
  for(i in 1:g){
    s <- q==i
    nobs <- sum(s); ne <- sum(outcome_num[s])
    
    dummystrat <- as.factor(rep("1", nobs))
    if(is.null(weights)){
      f = survfit(Srv[s,] ~ 1)
      pred[i] <- mean(predictions[s], na.rm=TRUE) #mean of predictions of each group
      ses_mean[i] <- sd(predictions[s], na.rm=TRUE)/sqrt(length(predictions[s])) #serror of predictions of each group
      
    }else{
      f = survfit(Srv[s,] ~ 1,weights =weights[s])
      pred[i] <- weighted.mean(predictions[s],w = weights[s], na.rm=TRUE) #mean of predictions of each group
      ses_mean[i] <- sqrt(wtd.var(predictions[s], w = weights[s],na.rm=TRUE))*sqrt(sum((weights[s]/sum(weights[s]))^2)) #serror of predictions of each group is the standard deviation divided by the effective n
      
    }
    
    #Note: if last u> last time, we take the estimates at the last time point:
    if(max(f$time)<u){
      print(paste0("Group",i," did not have any observation with survival time later than the end time point! The longest survival estimates were extended."))
      i.add = which.max(f$time)
      f$time = c(f$time,u)
      f$surv = c(f$surv,f$surv[i.add])
      f$std.err = c(f$std.err,f$std.err[i.add])
    }
    
    ##doesn't need conf.int since only need s.e.
    tt <- c(0, f$time)
    ss <- c(1, f$surv)
    se <- c(0, f$std.err)
    tm <- max((1:length(tt))[tt <= u+1e-6])
    km[i] <- ss[tm]
    std.err[i] <- se[tm]
    numobs[i]  <- nobs
    events[i]  <- ne
    n <- length(tt)
    if(u > tt[n]+1e-6 & ss[n]>0)
    {
      km[i] <- NA
      std.err[i] <- NA
    }
  }
  
  # Summarize all results:
  z <- cbind(x=pred, n=numobs, events=events, KM=km, 
             std.err=std.err)
  
  # Compute confidence intervals:
  # Note: This code was copied from the rms package. 
  # Essentially, it uses the fact that the margin of error is zcrit*std.err
  # the reason why we use the exp(zcrit*std.err) is because the std error is given for the log(st)
  
  ciupper <- function(surv, d) ifelse(surv==0, 0, pmin(1, surv*exp(d)))
  cilower <- function(surv, d) ifelse(surv==0, 0, surv*exp(-d))
  prop = km
  names(prop)=pred
  
  conf.int=0.95
  zcrit <- qnorm((conf.int+1)/2)
  low <- cilower(km, zcrit*std.err)
  hi  <- ciupper(km, zcrit*std.err)
  
  #print(f$lower[1])
  #print(f$surv[1]*exp(-qnorm((conf.int+1)/2)*f$std.err[1]))
  #stopifnot(round(f$lower[1],3)==round(f$surv[1]*exp(-qnorm((conf.int+1)/2)*f$std.err[1]),3)) # from the documentation in the survival package, f$std.err might be std error of logS or S, so we have this additional check
  
  
  # Plot confidence intervals:
  if(pl){
    if(log_plot){
      hi = ifelse(hi==1,0.999,hi)
    }
    if (surv_or_risk=="risk"){
      errbar(1-pred, 1-km, 1-hi, 1-low, add=TRUE)
    }else{
      errbar(pred, km, hi, low, add=TRUE)
    }
    
    
    #Add dots with predicted survival estimates and KM estimates
    if(surv_or_risk =="risk"){
      points(1-pred, 1-prop, pch=16,col=dots_col)
    }else{
      points(pred, prop, pch=16,col=dots_col)
    }
  }
  
  leg <- c(leg, "Grouped observations")
  lt <- c(lt, 0)
  col <- c(col, dots_col); lwd <- c(lwd, 1)     
  marks <- c(marks, 16)
  
  #Add probability distribution:
  if(surv_or_risk =="risk"){
    x<-1-predictions
  }else{
    x<-predictions
  }
  
  # Segments
  bins <- seq(lim[1], lim[2], length=101)
  x <- x[x >= lim[1] & x <= lim[2]]
  f <- table(cut(x, bins))
  j <- f > 0
  bins <- (bins[-101])[j]
  f <- f[j]
  
  # Add segments to plots:
  if(show_segments & pl){
    if(log_plot){ 
      #f <- lim[1] + .00015 * diff(lim) * f / max(f)
      #segments(bins, lim[1], bins, f) #This doesn't help because it seems like the probabilities are not concentrated below 0.1 when in fact they are
    }else{
      f <- lim[1] + .15 * diff(lim) * f / max(f)
      segments(bins, lim[1], bins, f,col =dots_col)
    }
  }
  
  # We return all of these quantities, to allow estiamtions of standard errors when pooling calibration plots across multiple NCC cohorts.
  return(list("logKMest_per_group"=log(km), "logKMse_per_group"=std.err,"Group_mean"=pred,"Serr_mean"=ses_mean))
}

# Pooled calibration plot 
#---------------------------------------------------
# @param cp_list: output list from the calib_plot function, with log of KM estimates, corresponding standard errors and mean of the risk predictions for each group in the plot.
# @param log_plot: TRUE if axes should be in log scale, FALSE otherwise
# @param lim: Vector of length 2, indicating the limits of the axes, e.g. c(0,1)
# @param surv_or_risk: String with two possible values: "surv" if survival probabilities should be shown, or "risk" if risk probabilities should be shown.
# @param u: number indicating at which timepoint the data should be truncated at.  It should be between 0 and max of outcome_FUP.
# @param dots_col: String indicating the color of the dots in the calibration plot
# @param file_name: string with the name of the file which will contain the  calibration plot
# @param main_plot: string with the title of the plot
calib_plot_pooled <- function(cp_list,log_plot,lim,surv_or_risk,u,dots_col="darkred",file_name=NULL,main_plot=""){
  
  if(!is.null(file_name)){
    pdf(file_name,width=12,height=8)
  }
  # Create a new plot:
  # Axes names:
  if(surv_or_risk =="risk"){
      ylab <- paste0("1-Kaplan-Meier survival ",u,' years') 
      xlab <- "Predicted risk probability"
    }else if (surv_or_risk =="surv"){
      ylab <- paste0("Kaplan-Meier survival ",u,' years') 
      xlab <- "Predicted survival probability"
  }
  
  if(log_plot){
      # Axes limits cannot contain 0
      if (lim[1]==0){
        print("Warning: Axes limits cannot contain 0 in log plots. Lower limit was replaced by 0.0001")
        lim[1] <-0.0001
      }
      plot(.5, .5, xlim=lim, ylim=lim, type="n", xlab=xlab, ylab=ylab,log ="xy",cex.lab=1.5, cex.axis=1.5,main = main_plot)
  }else{
      plot(.5, .5, xlim=lim, ylim=lim, type="n", xlab=xlab, ylab=ylab,cex.lab=1.5, cex.axis=1.5,main = main_plot)
  }
    abline(0, 1, lwd=6, col=gray(.85))
  
  # Extract points positions:
  km_sts <- exp(cp_list[["logKMest_per_group"]])
  km_sts[km_sts==0]<-0.000001 # This avoids NA values in the confidence intervals.
  qts <- log(-log(1-km_sts)) # Note that we pool variances based on complementary log log transformation, as suggested in https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-015-0048-4
  q_mean <- rowMeans(qts,na.rm=TRUE)
  km <- 1-exp(-exp(q_mean))
  m_groups <-cp_list[["Group_mean"]]
  m_pred <-rowMeans(m_groups)
  sd_pred <-cp_list[["Serr_mean"]]

  # Compute total variance, by taking into account both within and between variance.  
  # T_var <- U_m + (1+1/total_n)*B
  total_n <-ncol(km_sts)
  
  # Between variance:
  B <- apply(qts,1,stats::var)
  
  km_std <- cp_list[["logKMse_per_group"]]
  
  # We only have the variance of log KM estimates, so we use the delta method to write within variance as a function of var lox x.
  # Based on https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-015-0048-4
  U_m <-rowMeans((((km_std)^2*(km_sts)^2)/(log(1-km_sts)*(1-km_sts))^2),na.rm=TRUE)
  T_var <- U_m + (1+1/total_n)*B
  std.errs <-sqrt(T_var)
  
  # Add confidence intervals for the KM estimates:
  conf.int=0.95
  zcrit <- qnorm((conf.int+1)/2)
  
  low <- 1-exp(-exp(q_mean - zcrit*std.errs))
  hi  <- 1-exp(-exp(q_mean + zcrit*std.errs))
  

  if(log_plot){
      hi = ifelse(hi==1,0.999,hi)
    }
  if (surv_or_risk=="risk"){
      errbar(1-m_pred, 1-km, 1-hi, 1-low, add=TRUE)
    }else{
      errbar(m_pred, km, hi, low, add=TRUE)
   }

  #Add horizontal confidence intervals, for the variance regarding risk probabilities:
  
  B_h <- apply(m_groups,1,stats::var)
  U_m_h <- rowMeans(sd_pred^2)
  T_var_h <- U_m_h + (1+1/total_n)*B_h
  std.errs_h <-sqrt(T_var_h)
  low <- m_pred - zcrit*std.errs_h
  hi  <- m_pred + zcrit*std.errs_h
  
  if(log_plot){
    hi = ifelse(hi==1,0.999,hi)
  }
  
  for(j in c(1:length(low))){
    if (surv_or_risk=="risk"){
      arrows(x0=1-low[j], y0=1-km[j], x1=1-hi[j], y1=1-km[j], code=3, col="black", lwd=2,angle=90, length=0.05)
    }else{
      arrows(x0=0.1, y0=km[1], x1=0.2, y1=km[1], code=3, col="black", lwd=2,angle=90, length=0.5)
    }
  }

  #Add dots with predicted survival estimates and KM estimates
  if(surv_or_risk =="risk"){
    points(1-m_pred, 1-km, pch=16,col=dots_col)
  }else{
    points(m_pred, km, pch=16,col=dots_col)
  }
  if (!is.null(file_name)){
    dev.off()
  }                                                          
}

# Decision curve analysis:
#---------------------------------------------------
# Function to compute net benefit
# @param threshold: probability threshold
# @param outcome_num: numeric vector with status indicator, 0 =alive, 1 = dead. Follows the same conventions as parameter event in the Surv function from the survival package.
# @param outcome_FUP: numeric vector with follow up time of subjects
# @param risk_prob: numeric vector with model predictions
# @param tp: number indicating at which timepoint the data should be truncated at.  It should be between 0 and max of outcome_FUP.
# @param weights: if required output is weighted net benefit, weights is a numeric vector with weights of each observation. If required output is unweighted net benefit, weights should be NULL
computeNetBenefit <- function(threshold,outcome_num,outcome_FUP,risk_prob,tp,weights){
  #Truncate outcome at timepoint t:
  if (!is.null(tp)){
    truncated_survival <- truncateFUP(outcome_FUP,outcome_num,tp)
    outcome_FUP <- truncated_survival[["FUP"]]
    outcome_num <- truncated_survival[["Status"]]
  }
  
  total_n = length(outcome_num)
  
  # Compute the proportion of subjects with risk probability larger than the threshold:
  p_x1 =  weighted.mean(ifelse(risk_prob>threshold,1,0),weights)
  
  # Focus on high risk subjects:
  hr_pts = which(risk_prob>threshold)
  ev_tp_hr = outcome_FUP[hr_pts]
  ev_hr = outcome_num[hr_pts]
  weights_hr = weights[hr_pts]

  # Compute observed survival for high-risk subjects
  if(length(hr_pts)>1){
    expfit <- summary(survfit(Surv(ev_tp_hr, ev_hr) ~ 1,weights = weights_hr),times=tp,extend = TRUE)$surv
  }else{
    expfit = ifelse(ev_hr[hr_pts]==1,0,1)
  }
  
  # Compute true positives and false positives. Note True positives = trp* total_n and False positives = fp*total_n. However, these would be divided by total_n in the net benefit computation. So the total_n would cancel out.
  trp = (1-expfit)*p_x1 # trp = (1-expfit)*p_x1*total_n
  fp = expfit*p_x1 # fp = expfit*p_x1*total_n
  
  # Compute net benefit at the threshold 
  nb = trp-(fp)*(threshold/(1-threshold)) #nb = trp/total_n-(fp/total_n)*(threshold/(1-threshold))
  
  return(nb)
}

# Function to plot decision curve analysis
# @param ds: dataframe containing risk probabilities to be evaluated with a decision curve. Each column in this dataframe should correspond to a different risk prediction method
# @param vars: character vector with column names of columns present in dataframe ds, which represent risk probabilities of different methods
# @param col_vars: character vector with the color id for each variable present in the vars vector
# @param outcome_num: numeric vector with status indicator, 0 =alive, 1 = dead. Follows the same conventions as parameter event in the Surv function from the survival package.
# @param outcome_FUP: numeric vector with follow up time of subjects
# @param ylim: vector of length two, with lower and upper limit for the y axis in the plot
# @param xlim: vector of length two, with lower and upper limit for the x axis in the plot
# @param tp: number indicating at which timepoint the data should be truncated at.  It should be between 0 and max of outcome_FUP.
# @param weights: if required output is weighted decision curve, weights is a numeric vector with weights of each observation. If required output is unweighted decision curve, weights should be NULL
# @param var_labels: vector with the labels that should be used for each variable. Vector names should correspond to vars.
# @param plot: boolean indicating whether decision curve should be displayed on a plot (TRUE) or not (FALSE)
dca_survival = function(ds,vars,col_vars,outcome_num,outcome_FUP,ylim,xlim,tp,weights,var_labels=NULL,plot=TRUE){
  
  #Store all net benefits
  list_nbs <-vector("list",1+length(vars))
  names(list_nbs) <-c("NB_all",paste0("NB_",vars))
  
  nb_th = seq(0,0.99,0.01)
  #Net benefit of treating all:
  list_nbs[["NB_all"]]<-nb_all<-sapply(nb_th,computeNetBenefit,outcome_num,outcome_FUP,rep(1,length(outcome_num)),tp,weights)
  if(plot){
    plot(nb_th*100,nb_all,type="l",ylim=ylim,xlim=xlim,ylab = "Net benefit",xlab="Threshold probability (%)",col="black",lwd=2,cex.lab=1.5, cex.sub=1.5,cex.axis=1.2)
  }

  
  #Net benefit of methods of interest:
  for (var in vars){
    list_nbs[[paste0("NB_",var)]]<-nb_var<-as.numeric(sapply(nb_th,computeNetBenefit,outcome_num,outcome_FUP,ds[,var],tp,weights))
    if(plot){
      lines(nb_th*100,nb_var,lwd=2,col=col_vars[var])
    }
    
  }
  
  #Net benefit of treating none:
  if (plot){
    lines(nb_th*100,rep(0,length(nb_th)),col="darkgray",lwd=2)
  }
  
  
  #Add a legend:
  if (plot){
    if(is.null(var_labels)){
      legend("topright",legend=c("Treat all","Treat none",vars),lty = 1,lwd=2, col = c("black","darkgray",col_vars),bty="n",cex=1.2)
    }else{
      legend("topright",legend=c("Treat all","Treat none",var_labels),lty = 1,lwd=2, col = c("black","darkgray",col_vars),bty="n",cex=1.2)
    }
  }
  return(list_nbs)
}

# Function to plot decision curve analysis, with confidence intervals for the realizations of the NCC datasets
# @param list_nbs: list with net benefits to be plotted
# @param vars: character vector with column names of columns present in dataframe ds, which represent risk probabilities of different methods
# @param col_vars: character vector with the color id for each variable present in the vars vector
# @param ylim: vector of length two, with lower and upper limit for the y axis in the plot
# @param xlim: vector of length two, with lower and upper limit for the x axis in the plot
# @param var_labels: vector with the labels that should be used for each variable. Vector names should correspond to vars.
# @param file_name: string with the name of the file which will contain the decision curve plot
dca_survival_plot_ci = function(list_nbs,vars,col_vars,ylim,xlim,var_labels=NULL,file_name=NULL){
  
  if (!is.null(file_name)){
    pdf(file_name,width=12,height=8)
  }
  nb_th = seq(0,0.99,0.01)
  
  #Net benefit of treating all:
  nb_all<-rowMeans(list_nbs[["NB_all"]],na.rm=TRUE)
  plot(nb_th*100,nb_all,type="l",ylim=ylim,xlim=xlim,ylab = "Net benefit",xlab="Threshold probability (%)",col="black",lwd=2,cex.lab=1.5, cex.sub=1.5,cex.axis=1.2)
  
  nb_all_qnts <-apply(list_nbs[["NB_all"]],1,function(x) c(quantile (x,0.025),quantile (x,0.975)))
  #lines(nb_th*100,nb_all_qnts[1,],lwd=2,lty=2)
  #lines(nb_th*100,nb_all_qnts[2,],lwd=2,lty=2)
  polygon(c(nb_th*100,rev(nb_th*100)),c(nb_all_qnts[1,],rev(nb_all_qnts[2,])),col= adjustcolor("black", alpha.f = 0.2), border = FALSE)
  
  
  #Net benefit of methods of interest:
  for (var in vars){
    nb_var <- rowMeans(list_nbs[[paste0("NB_",var)]],na.rm=FALSE) 
    
    nb_vars_qnts <-apply(list_nbs[[paste0("NB_",var)]],1,function(x) c(quantile (x,0.025,na.rm=TRUE),quantile (x,0.975,na.rm=TRUE)))
    nb_vars_qnts_nas = is.na(nb_var)
    polygon(c(nb_th[!nb_vars_qnts_nas]*100,rev(nb_th[!nb_vars_qnts_nas]*100)),c(nb_vars_qnts[1,!nb_vars_qnts_nas],rev(nb_vars_qnts[2,!nb_vars_qnts_nas])),col= adjustcolor(col_vars[var], alpha.f = 0.2), border = FALSE)
   
    lines(nb_th*100,nb_var,lwd=2,col=col_vars[var])
    #lines(nb_th*100,nb_vars_qnts[1,],lwd=2,lty=2,col=col_vars[var])
    #lines(nb_th*100,nb_vars_qnts[2,],lwd=2,lty=2,col=col_vars[var])
  }
  
  #Net benefit of treating none:
  lines(nb_th*100,rep(0,length(nb_th)),col="darkgray",lwd=2)
  
  #Add a legend:

  if(is.null(var_labels)){
      legend("topright",legend=c("Treat all","Treat none",vars),pch = 16, col = c("black","darkgray",col_vars),bty="n")
  }else{
      legend("topright",legend=c("Treat all","Treat none",var_labels),pch = 16, col = c("black","darkgray",col_vars[var_labels]),bty="n")
  }
  
  if (!is.null(file_name)){
    dev.off()
  }
}

