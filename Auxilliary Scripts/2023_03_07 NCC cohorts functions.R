#---------------------------------------------------
# Description:Auxilliary functions for the derivation of nested case-control cohorts
# Author: B. Rentroia Pacheco
#---------------------------------------------------

# Create a Nested Case-control (NCC) cohort for prediction model development/validation
# For these scenarios, cases should not be selected as controls for other cases. Controls should not be selected twice.
#---------------------------------------------------
# @param outcome_FUP_name: name of the column with follow up time of subjects
# @param outcome_num_name: name of the column with status indicator, 0 =alive, 1 = dead. Follows the same conventions as parameter event in the Surv function from the survival package.
# @param ds: full cohort, with columns containing outcome indicator and outcome follow up, matching variables and any variables that should be kept for the final NCC dataframe
# @param include_vars: vector with column names that should be kept in the final NCC dataframe
# @param match_vars: vector with column names of variables that will be used for matching controls to cases. Note that ccwc only deals with categorical variables for matching.
# @param ncontrols: number of controls per case
# @param pt_id: name of the column that contains the id of the patients.
# @param seed: number of the random seed
# @param tp: number indicating at which timepoint at which prevalence should be computed. It should be between 0 and max of outcome_FUP.
NCCforPredModel <-function (outcome_FUP_name,outcome_num_name,ds,include_vars,match_vars,ncontrols,pt_id,wgh_method,seed,tp=NULL){
  
  fullcohort = ds
  
  # Initial NCC dataset:
  ncc_cohort = ccwc2(exit=eval(sym(outcome_FUP_name)), fail=eval(sym(outcome_num_name)),origin =0, include = include_vars,match = match_vars, data=ds, controls=ncontrols, silent=TRUE,seed = seed)
  
  extend_ncc = TRUE
  while (extend_ncc){
    # Cases should not be selected as controls for other cases. Controls should not be selected twice.
    cases_ncc = ncc_cohort[which(ncc_cohort[,"Fail"]==1),pt_id]
    controls_ncc = ncc_cohort[which(ncc_cohort[,"Fail"]==0),pt_id]
    
    #Remove cases that were considered controls as potential controls:
    cases_that_were_considered_controls = intersect(cases_ncc,controls_ncc)
    if (length(cases_that_were_considered_controls)>0){
      ncc_cohort = ncc_cohort %>% filter(!(!!rlang::sym(pt_id)%in%cases_that_were_considered_controls & Fail ==0)) %>%as.data.frame()
    }
    
    # Keep only one random record for controls that were selected multiple times:
    controls_repeated = names(table(controls_ncc)[table(controls_ncc)>1])
    if (length(controls_repeated)>0){
      ncc_cohort = ncc_cohort %>% group_by(!!rlang::sym(pt_id)) %>% sample_n(1) %>%ungroup() %>% as.data.frame()
    }
    
    incomplete_sets = names(table(ncc_cohort[,"Set"])[table(ncc_cohort[,"Set"])<2])
    
    if(length(incomplete_sets)==0){
      extend_ncc <- FALSE
    }else if (length(incomplete_sets)>0){
      # Repeat assignment of controls to cases:
      cases_and_controls_with_pair = ncc_cohort[!ncc_cohort[,"Set"]%in%incomplete_sets,pt_id]
      ds = ds %>% filter(!(!!rlang::sym(pt_id)%in% cases_and_controls_with_pair)) %>%as.data.frame()
      ncc_cohort_2 = ccwc2(exit=eval(sym(outcome_FUP_name)), fail=eval(sym(outcome_num_name)),origin =0, data=ds, include=include_vars, match = match_vars ,controls=ncontrols, silent=TRUE,seed = seed)
      ncc_cohort = ncc_cohort %>% filter(!(Set%in% incomplete_sets)) %>%as.data.frame()
      
      # Merge NCC cohort:
      ncc_cohort$Set = as.numeric(factor(ncc_cohort$Set))
      ncc_cohort_2$Set = as.numeric(factor(ncc_cohort_2$Set))+max(ncc_cohort$Set)
      ncc_cohort = rbind(ncc_cohort,ncc_cohort_2)
    }
    
  }
  
  # Add sampling weights:
  ncc_cohort<-computeSamplingWeights(fullcohort[,outcome_FUP_name],fullcohort[,outcome_num_name],fullcohort,ncc_cohort ,pt_id,ncontrols,wgh_method,match_vars,tm=tp)
  return(ncc_cohort)
}

# Creates a typical NCC cohort.
# The function ccwc from the EPI package was modified to:
# 1. allow more flexibility in the variables provided for the parameters "include" and "match"
# 2. set seed for reproducibility
#---------------------------------------------------
# @param entry: Time of entry to follow-up
# @param exit: Time of exit from follow-up Status
# @param fail:  Status on exit (1=Fail, 0=Censored)
# @param origin: Origin of analysis time scale
# @param controls: The number of controls to be selected for each case
# @param match: Vector of categorical variables on which to match cases and controls
# @param include: Vector of other variables to be carried across into the case-control study
# @param data: Data frame in which to look for input variables
# @param silent: If FALSE, echos a . to the screen for each case-control set created; otherwise produces no output.
# @param seed: number of the random seed
ccwc2 <-
  function(entry=0, exit, fail, origin=0, controls=1, match=NULL,
           include=NULL, data=NULL, silent=FALSE,seed = 123){
    
    set.seed(seed) # for reproducibility
    # Check arguments
    entry <- eval(substitute(entry), data)
    exit <- eval(substitute(exit), data)
    fail <- eval(substitute(fail), data)
    origin <- eval(substitute(origin), data)
    
    n <- length(fail)
    if (length(exit)!=n)
      stop("All vectors must have same length")
    if (length(entry)!=1 && length(entry)!=n)
      stop("All vectors must have same length")
    if (length(origin)==1) {
      origin <- rep(origin, n)
    }
    else {
      if (length(origin)!=n)
        stop("All vectors must have same length")
    }
    # Transform times to correct scale
    t.entry <- as.numeric(entry - origin)
    t.exit <- as.numeric(exit - origin)
    
    # match= argument
    #marg <- substitute(match)
    #if (mode(marg)=="name") {
    #  match <- list(eval(marg, data))
    #  names(match) <- as.character(marg)
    #}
    #else if (mode(marg)=="call" && marg[[1]]=="list") {
    #  mnames <- names(marg)
    #  nm <- length(marg)
    #  if (nm>1) {
    #    if (!is.null(mnames)) {
    #      for (i in 2:nm) {
    #        if (mode(marg[[i]])=="name")
    #          mnames[i] <- as.character(marg[[i]])
    #        else
    #          stop("illegal argument (match)")
    #      }
    #    }
    #    else {
    #      for (i in 2:nm) {
    #        if (mode(marg[[i]])=="name")
    #          mnames[i] <- as.character(marg[[i]])
    #        else
    #          stop("illegal argument (match)")
    #      }
    #      mnames[1] <= ""
    #    }
    #  }
    #  names(marg) <- mnames
    #  match <- eval(marg, data)
    #}
    #else {
    #  stop("illegal argument (match)")
    #}
    
   if (length(match)==1){
      match_list = list(data[,match])
      names(match_list) = match
      match=match_list
    } else if (length(match)>1){
     match = as.list(data[,match])
   }
    
    m <- length(match)
    mnames <- names(match)
    if (m>0) {
      for (i in 1:m) {
        if (length(match[[i]])!=n) {
          stop("incorrect length for matching variable")
        }
      }
    }
    # include= argument
    #iarg <- substitute(include)
    #if (mode(iarg)=="name") {
    #  include <- list(eval(iarg, data))
    #  names(include) <- as.character(iarg)
    #}
    #else if (mode(iarg)=="call" && iarg[[1]]=="list") {
    #  ni <- length(iarg)
    #  inames <- names(iarg)
    #  if (ni>1) {
    #    if (!is.null(inames)) {
    #      for (i in 2:ni) {
    #        if (mode(iarg[[i]])=="name")
    #          inames[i] <- as.character(iarg[[i]])
    #        else
    #          stop("illegal argument (include)")
    #      }
    #    }
    #    else {
    #      for (i in 2:ni) {
    #        if (mode(iarg[[i]])=="name")
    #          inames[i] <- as.character(iarg[[i]])
    #        else
    #          stop("illegal argument (include)")
    #      }
    #      inames[1] <= ""
    #    }
    #  }
    #  names(iarg) <- inames
    #  include <- eval(iarg, data)
    #}
    #else {
    #  stop("illegal argument (include)")
    #}
    
   if (length(include)==1){
      include = list(data[,include])
    }else if(length(include)>1){
      include = as.list(data[,include])
    }
    ni <- length(include)
    inames <- names(include)
    if (ni>0) {
      for (i in 1:ni) {
        if (length(include[[i]])!=n) {
          stop("incorrect length for included variable")
        }
      }
    }
    # create group codes using matching variables
    grp <- rep(1,n)
    pd <- 1
    if (m>0) {
      for (im in 1:m) {
        v <- match[[im]]
        if (length(v)!=n)
          stop("All vectors must have same length")
        if (!is.factor(v))
          v <- factor(v)
        grp <- grp + pd*(as.numeric(v) - 1)
        pd <- pd*length(levels(v))
      }
    }
    # Create vectors long enough to hold results
    nn <- (1+controls)*sum(fail!=0)
    pr <- numeric(nn)
    sr <- numeric(nn)
    tr <- vector("numeric", nn)
    fr <- numeric(nn)
    nn <- 0
    # Sample each group
    if (!silent) {
      cat("\nSampling risk sets: ")
    }
    set <- 0
    nomatch <- 0
    incomplete <- 0
    ties <- FALSE
    fg <- unique(grp[fail!=0])
    for (g in fg) {
      # Failure times
      ft <- unique( t.exit[(grp==g) & (fail!=0)] )
      # Create case-control sets
      for (tf in ft) {
        if (!silent) {
          cat(".")
        }
        set <- set+1
        case <- (grp==g) & (t.exit==tf) & (fail!=0)
        ncase <- sum(case)
        if (ncase>1)
          ties <- TRUE
        noncase <- (grp==g) & (t.entry<=tf) &
          (t.exit>=tf) & !case
        ncont <- controls*ncase
        if (ncont>sum(noncase)) {
          ncont <- sum(noncase)
          if (ncont>0) incomplete <- incomplete + 1
        }
        if (ncont>0) {
          newnn <- nn+ncase+ncont
          sr[(nn+1):newnn] <- set
          tr[(nn+1):newnn] <- tf
          fr[(nn+1):(nn+ncase)] <- 1
          fr[(nn+ncase+1):newnn] <- 0
          pr[(nn+1):(nn+ncase)] <- (1:n)[case]
          ## Work around bad behaviour of sample for vectors of length 1
          noncase.id <- (1:n)[noncase]
          pr[(nn+ncase+1):(newnn)] <- noncase.id[sample.int(length(noncase.id),
                                                            size=ncont)]
          nn <- newnn
        }
        else {
          nomatch <- nomatch + ncase
        }
      }
    }
    if (!silent) {
      cat("\n")
    }
    res <- vector("list", 4+m+ni)
    if (nn>0) {
      res[[1]] <- sr[1:nn]
      res[[2]] <- map <- pr[1:nn]
      res[[3]] <- tr[1:nn] + origin[map]
      res[[4]] <- fr[1:nn]
    }
    if (m>0) {
      for (i in 1:m) {
        res[[4+i]] <- match[[i]][map]
      }
    }
    if (ni>0) {
      for (i in 1:ni) {
        res[[4+m+i]] <- include[[i]][map]
      }
    }
    names(res) <- c("Set", "Map", "Time", "Fail", mnames, inames)
    if (incomplete>0)
      warning(paste(incomplete, "case-control sets are incomplete"))
    if (nomatch>0)
      warning(paste(nomatch, "cases could not be matched"))
    if (ties)
      warning("there were tied failure times")
    data.frame(res)
  }

# Compute Sampling weights
# The multipleNCC package estimates sampling probabilities using different methods. 
# They require a vector containing sampling and status information, for the entire cohort.
#---------------------------------------------------
# @param outcome_FUP: numeric vector with follow up time of subjects
# @param outcome_num: numeric vector with status indicator, 0 =alive, 1 = dead. Follows the same conventions as parameter event in the Surv function from the survival package.
# @param ds: full cohort, with additional variables
# @param ncc_ds: dataframe containing the nested case cohort
# @param pt_id: name of the column that contains the id of the patients.
# @param ncontrols: number of controls per case
# @param w_method: character vector describing the weighting method for the calculation of sampling weights
# @param match_vars: vector with column names of variables that will be used for matching controls to cases. Note that ccwc only deals with categorical variables for matching.
# @param tm: number indicating at which timepoint at which prevalence should be computed.  It should be between 0 and max of outcome_FUP.
computeSamplingWeights <-function(outcome_FUP,outcome_num,ds,ncc_ds,pt_id,ncontrols,w_method,match_vars,tm = NULL){
  samplestat <- outcome_num *2 # With this, we have controls = 0, cases (events) =2
  samplestat <- ifelse(ds[,pt_id]%in%ncc_ds[which(ncc_ds[,"Fail"]==0),pt_id],1,samplestat)
  
  if(!is.null(match_vars)){
    m.vars <- cbind(ds[,match_vars])
    m.int <- c(0,0)
  }else{
    m.vars <- 0
    m.int <-0
  }
  
  # Methods provided by the multipleNCC package:
  if(w_method =="KM"){
    p_samp <- KMprob(outcome_FUP, samplestat, ncontrols, left.time = 0, match.var = m.vars, match.int = m.int)
  }else if (w_method =="GAMprob"){
    p_samp <- GAMprob(outcome_FUP, samplestat, left.time = 0, match.var = m.vars, match.int = m.int)
  }else if (w_method =="GLMprob"){
    p_samp <- GLMprob(outcome_FUP, samplestat, left.time = 0, match.var = m.vars, match.int = m.int)
  }else if (w_method =="Chenprob"){
    # The multipleNCC package provides the possibility of computing Chen weights, but only in situations without matching, so we are not going to include these weights in the experiments.  
  }
  
  if(!w_method%in%c("Prevalence","Prevalence_KM")){
    sampling_prob_df<-data.frame(pt_id = ds[,pt_id],samp_prob = p_samp,samp_weight = 1/p_samp) 
    
    final_ds <- merge(ncc_ds,sampling_prob_df,by.x = pt_id,by.y="pt_id")
# Naive methods that use only prevalence of the event in the cohort:
    
  }else{
    final_ds <-ncc_ds
    final_ds$samp_prob <- 1
    prevalence <- 1-summary(survfit(Surv(ds[,outcome_FUP_name],ds[,outcome_num_name])~1),times=tm)$surv
    nr_controls <- sum(final_ds$Fail==0)
    
    if (w_method =="Prevalence"){
      final_ds$samp_prob[which(final_ds$Fail==0)] <- nr_controls/(sum(ds[,outcome_num_name]==0)) # this is just the fraction of controls that are selected
    }else if (w_method =="Prevalence_KM"){
      final_ds$samp_prob[which(final_ds$Fail==0)]  <- 1- (1-1/sum(ds[,outcome_num_name]))^nr_controls # computed as the complement of the probability of never being sampled.
    }
    final_ds$samp_weight<-1/final_ds$samp_prob
    
  }
  
  # Rescale weights
  tot_controls <- sum(ds[,outcome_num_name]==0)
  final_ds$samp_weight[which(final_ds$Fail==0)] <- final_ds$samp_weight[which(final_ds$Fail==0)]*(tot_controls/sum(final_ds$samp_weight[which(final_ds$Fail==0)]))
  
  
  return(final_ds)
}


# Wrapper function to evaluate model predictions in 1 NCC cohort
#---------------------------------------------------
# @param outcome_FUP: numeric vector with follow up time of subjects
# @param outcome_num: numeric vector with status indicator, 0 =alive, 1 = dead. Follows the same conventions as parameter event in the Surv function from the survival package.
# @param ds: full cohort, with additional variables
# @param include_vars: vector with column names that should be kept in the final NCC dataframe
# @param match_vars: vector with column names of variables that will be used for matching controls to cases. Note that ccwc only deals with categorical variables for matching.
# @param ncontrols: number of controls per case
# @param pt_id: name of the column that contains the id of the patients.
# @param wgh_method: character vector describing the weighting method for the calculation of sampling weights
# @param seed: number of the random seed
# @param pt_id: name of the column that contains the model risk predictions
# @param tp: number indicating at which timepoint the data should be truncated at.  It should be between 0 and max of outcome_FUP.
# @param nrep: number of derivations of NCC cohorts
# @param ylim_dca: vector of length two, with lower and upper limit for the y axis in the decision curve plot
# @param xlim_dca: vector of length two, with lower and upper limit for the x axis in the decision curve plot
# @param ylim_dca_un_fact: factor that will be used to multiply the provided limits of the y axis in the unweighted decision curve plot
# @param file_name_dca: string with the name of the file which will contain the decision curve plot
# @param cp_prob_range: Vector of length 2, indicating the limits of the axes, e.g. c(0,1)
# @param cp_nr_pt_gr: Number of groups that will be displayed in the plot
# @param ylim_cp_un_fact: factor that will be used to multiply the provided limits of the y axis in the unweighted calibration plot
# @param file_name_cp: string with the name of the file which will contain the calibration plot
# @param cp_fullcohort: list that is the output of the calib_plot for the model applied to the full cohort. This is used to add the calibration points obtained in the full cohort to the pooled calibration plot.
# @param dc_fullcohort: list that is the output of the dca_survival for the model applied to the full cohort. This is used to add the decision curve obtained in the full cohort to the pooled decision curve plot.
evaluateModelInNCC  <-function (outcome_FUP_name,outcome_num_name,ds,include_vars,match_vars,ncontrols,pt_id,wgh_method,seed,predictions_name,tp,nrep,ylim_dca,xlim_dca,ylim_dca_un_fact,file_name_dca,cp_prob_range,cp_nr_pt_gr,ylim_cp_un_fact,file_name_cp,cp_fullcohort,dc_fullcohort){
  
  # Data structure to store performance metrics:
  performance_metrics_ncc <- vector("list", 14)
  names(performance_metrics_ncc)<- c("un_Cind","un_OE","un_CSlope","un_SE","un_SP","un_PPV","un_NPV","wgh_Cind","wgh_OE","wgh_CSlope","wgh_SE","wgh_SP","wgh_PPV","wgh_NPV")
  performance_metrics_ncc <- lapply(performance_metrics_ncc,function(x) vector("numeric",nrep))

  # And for net benefits: 
  nb_list_wgh<- list(as.data.frame(matrix(NA, ncol = nrep,nrow =100)),as.data.frame(matrix(NA, ncol = nrep,nrow =100)))
  names(nb_list_wgh) <-c("NB_all",paste0("NB_",predictions_name))
  nb_list_un <- nb_list_wgh
  
  # And for calibration plots:
  cp_list_wgh<- list(as.data.frame(matrix(NA, ncol = nrep,nrow =cp_nr_pt_gr)),as.data.frame(matrix(NA, ncol = nrep,nrow =cp_nr_pt_gr)),as.data.frame(matrix(NA, ncol = nrep,nrow =cp_nr_pt_gr)),as.data.frame(matrix(NA, ncol = nrep,nrow =cp_nr_pt_gr)))
  names(cp_list_wgh) <-c("logKMest_per_group","logKMse_per_group","Group_mean","Serr_mean")
  cp_list_un <- cp_list_wgh
  
  
  col_vars <- "darkred"
  names(col_vars) <-predictions_name
  
  # Always include matching variables in the ncc cohort:
  if (!is.null(match_vars)){
    include_vars <-c(include_vars,match_vars)
  }
  
  # Repeat nrep times:
  for (i in c(1:nrep)){
    
    # Generate a NCC cohort, and corresponding sampling weights:
    ncc_cohort <- NCCforPredModel(outcome_FUP_name,outcome_num_name,ds=ds,include = include_vars ,match_vars = match_vars,ncontrols,pt_id,wgh_method, seed=seed+i,tp=tp)
     
    # Compute weighted and unweighted metrics:
    for (wgh in c("un","wgh")){
      if(wgh=="un"){
        weight_scheme <- NULL
        weight_scheme_dca <-rep(1,nrow(ncc_cohort))
      }else if(wgh=="wgh"){
        weight_scheme <- weight_scheme_dca <- ncc_cohort$samp_weight
        
      }
      
      performance_metrics_ncc[[paste0(wgh,"_Cind")]][i] <- computeCindex(ncc_cohort[,outcome_FUP_name],ncc_cohort[,"Fail"],ncc_cohort[,predictions_name],tp,weights = weight_scheme)
      performance_metrics_ncc[[paste0(wgh,"_OE")]][i]  <- computeOEratio(ncc_cohort[,outcome_FUP_name],ncc_cohort[,"Fail"],ncc_cohort[,predictions_name],tp,weights = weight_scheme)
      performance_metrics_ncc[[paste0(wgh,"_CSlope")]][i]  <- calibrationSlope(ncc_cohort[,outcome_FUP_name],ncc_cohort[,"Fail"],ncc_cohort[,predictions_name],tp,weights = weight_scheme)
      
      cutoff_perf_metrics <-computeDiscMetrics(ncc_cohort[,outcome_FUP_name],ncc_cohort[,"Fail"],ncc_cohort[,predictions_name],0.03,tp,weights=weight_scheme)
      
      for (mt in c(1:4)){
        performance_metrics_ncc[[paste0(wgh,c("_SE","_SP","_PPV","_NPV"))[mt]]][i]<- cutoff_perf_metrics[mt]
      }
      
      # Decision curves:
      nb_list_ncc <- dca_survival(ncc_cohort,predictions_name,col_vars ,ncc_cohort[,"Fail"],ncc_cohort[,outcome_FUP_name],ylim_dca,xlim_dca,tp,weight_scheme_dca,var_labels=NULL,plot=FALSE)
      if(wgh =="un"){
        nb_list_un[["NB_all"]][,i] <-nb_list_ncc[["NB_all"]]
        nb_list_un[[paste0("NB_",predictions_name)]][,i] <-nb_list_ncc[[paste0("NB_",predictions_name)]]
      }else if (wgh =="wgh"){
        nb_list_wgh[["NB_all"]][,i] <-nb_list_ncc[["NB_all"]]
        nb_list_wgh[[paste0("NB_",predictions_name)]][,i] <-nb_list_ncc[[paste0("NB_",predictions_name)]]
      }
      
      # Calibration plot:
      cp_col_dots <- ifelse(wgh=="un","#F8766D","#00BFC4")
      cp_list_ncc <-  calib_plot(ncc_cohort[,predictions_name],ncc_cohort[,"Fail"],ncc_cohort[,outcome_FUP_name],tp,floor(nrow(ncc_cohort)/cp_nr_pt_gr),lim=cp_prob_range,pl=FALSE,weights=weight_scheme,log_plot=FALSE,main_plot="",surv_or_risk = "risk",new_plot = FALSE,dots_col = cp_col_dots,show_segments=FALSE)      
      
      for(cln_cp in names(cp_list_ncc)){
        if(wgh =="un"){
          cp_list_un[[cln_cp]][,i] <- as.numeric(cp_list_ncc[[cln_cp]])
        }else if (wgh =="wgh"){
          cp_list_wgh[[cln_cp]][,i] <- as.numeric(cp_list_ncc[[cln_cp]])
        }
      }
      
    }
  }
  
  # Summarize all performance metrics, using mean+2.5 and 97.5% quantiles:
  summary_performance_ncc <- lapply(performance_metrics_ncc, function(x) paste0(round(mean(x),2),"(",round(quantile(x,0.025),2),"-",round(quantile(x,0.975),2),")"))
  summary_performance_ncc <- unlist(summary_performance_ncc)
  
  # Summarize decision curve analysis:
  col_vars_dca_plot <- "#00BFC4"
  names(col_vars_dca_plot) <-predictions_name
  dca_survival_plot_ci(nb_list_wgh,predictions_name,col_vars_dca_plot,ylim_dca,xlim_dca,var_labels=NULL,paste0(file_name_dca,"_wgh.pdf"))
  col_vars_dca_plot <- "#F8766D"
  names(col_vars_dca_plot) <-predictions_name
  dca_survival_plot_ci(nb_list_un,predictions_name,col_vars_dca_plot,ylim_dca*ylim_dca_un_fact,xlim_dca,var_labels=NULL,paste0(file_name_dca,"_unwgh.pdf"))
  
  # Summarize calibration plots:
  calib_plot_pooled(cp_list_wgh,log_plot=FALSE,lim=cp_prob_range,"risk",tp,dots_col = "#00BFC4",file_name=paste0(file_name_cp,"_wgh.pdf"))
  calib_plot_pooled(cp_list_un,log_plot=FALSE,lim=cp_prob_range*ylim_cp_un_fact,"risk",tp,dots_col = "#F8766D",file_name=paste0(file_name_cp,"_un.pdf"))
  
  # Combined plot:
  pdf(paste0(file_name_cp,"_allplots.pdf"),width=13,height=6)
  par(mfrow=c(1,2))
  calib_plot_pooled(cp_list_un,log_plot=FALSE,lim=cp_prob_range*ylim_cp_un_fact,"risk",tp,dots_col="#F8766D",file_name = NULL,main_plot = "Unweighted")
  # Add full cohort estimates:
  points(1-cp_fullcohort$Group_mean,1-exp(cp_fullcohort$logKMest_per_group),col = "blue",pch=16)
  legend("topright",legend=c("NCC","Full cohort"),pch=16, col = c("#F8766D","blue"),bty="n",cex=1.2)
  
  calib_plot_pooled(cp_list_wgh,log_plot=FALSE,lim=cp_prob_range,"risk",tp,dots_col="#00BFC4",file_name = NULL,main_plot = "Weighted")
  # Add full cohort estimates:
  points(1-cp_fullcohort$Group_mean,1-exp(cp_fullcohort$logKMest_per_group),col = "blue",pch=16)
  legend("topright",legend=c("NCC","Full cohort"),pch=16, col = c("#00BFC4","blue"),bty="n",cex=1.2)
  
  dev.off ()
  
  pdf(paste0(file_name_dca,"_allplots.pdf"),width=13,height=6)
  par(mfrow=c(1,2))
  # Add full cohort in the legend:
  col_vars <-c("#F8766D","#F8766D","blue")
  names(col_vars)<-c(predictions_name,"NCC","Full cohort")
  var_labels_dca_ci <- c("NCC","Full cohort")
  dca_survival_plot_ci(nb_list_un,predictions_name,col_vars,ylim_dca*ylim_dca_un_fact,xlim_dca,var_labels=var_labels_dca_ci,NULL)
  # Add full cohort estimates:
  lines(seq(0,0.99,0.01)*100,dc_fullcohort[[paste0("NB_",predictions_name)]],lwd=2,col="blue",lty=2)
  col_vars[c(predictions_name,"NCC")] <- rep("#00BFC4",2)
  dca_survival_plot_ci(nb_list_wgh,predictions_name,col_vars,ylim_dca,xlim_dca,var_labels=var_labels_dca_ci,NULL)
  # Add full cohort estimates:
  lines(seq(0,0.99,0.01)*100,dc_fullcohort[[paste0("NB_",predictions_name)]],lwd=2,col="blue",lty=2)
  
  dev.off()
  
  return(list("Summary_performance"=summary_performance_ncc,"Performance_metrics_all" = performance_metrics_ncc))  
}

# TESTS:
# ncc_cohort = NCCforPredModel(outcome_num_name,outcome_FUP_name,ds=ERGOcohort,include = include_vars ,match_vars = NULL,1,"ergoid",seed=123,wgh_method="KM")
# stopifnot(all(table(ncc_cohort$Fail)==163)) #We should have 163 cases and 163 controls
# ncc_cohort_2 = NCCforPredModel(outcome_num_name,outcome_FUP_name,ds=ERGOcohort,include = include_vars ,match_vars = NULL,2,"ergoid",seed=123,wgh_method="KM")
# stopifnot(all(table(ncc_cohort_2$Fail)==c(303,163))) We should have 163 cases and 163*2 controls
# stopifnot(all(ncc_cohort$TimeToEvent_max10-ncc_cohort$Time>=0)) #FU of controls should be always higher than cases
# stopifnot(all(duplicated(ncc_cohort$ergoid)==FALSE)) #No duplicates
# ncc_cohort = ccwc2(exit=eval(sym(outcome_FUP_name)), fail=eval(sym(outcome_num_name)),origin =0, include = c("ergoid","startage") ,data=ERGOcohort, controls=1, silent=TRUE,seed = 123,wgh_method="KM")
# stopifnot(any(duplicated(ncc_cohort$ergoid)==FALSE))
# Test matching variables:
# ncc_cohort = NCCforPredModel(outcome_num_name,outcome_FUP_name,ds=ERGOcohort,include = include_vars ,match_vars = "age_categories",1,"ergoid",seed=123,wgh_method="KM")
# stopifnot(all(names(table(table(ncc_cohort$Set[ncc_cohort$age_categories==">60"])))==2))
# Note: I checked sampling probabilities for a test example:
# ncc_cohort <- NCCforPredModel(outcome_FUP_name,outcome_num_name,ds=ds[1:200,],include = include_vars ,match_vars = "age_categories",1,"ergoid","KM", seed=123+i)
# Controls with follow up time =10 and <60, should have the following probability: 1-(115/116)*(116/117)*(118/119)*(120/121)
# stopifnot(all(ncc_cohort$samp_prob[which(ncc_cohort$Fail==0 & ncc_cohort$TimeToEvent_max10==10 & ncc_cohort$age_categories=="<60")] == (1-(115/116)*(116/117)*(118/119)*(120/121))))
# stopifnot(round(ncc_cohort$samp_prob[which(ncc_cohort$Fail==0 & ncc_cohort$TimeToEvent_max10==10 & ncc_cohort$age_categories==">60")],5)==round(1/69,5))

# Tests: weight computation:
# stopifnot(all(ncc_cohort$samp_weight[ncc_cohort$Fail==1]==1))
