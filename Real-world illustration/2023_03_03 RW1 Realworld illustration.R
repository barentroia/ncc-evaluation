#---------------------------------------------------
# Aim: Reproduce results of the publication: 
# Weighted metrics are required when evaluating the performance of prediction models in nested case-control cohorts
# Author: B. Rentroia Pacheco
# Input: Preprocessed datasets from the Lakeman 2020 publication
# Output: Results in the Rentroia-Pacheco 2023 publication.
#---------------------------------------------------

#---------------------------------------------------
# 0. Setup
#---------------------------------------------------
# Libraries:
library(survC1) # For computation of Uno's concordance metric
library(forcats)
library(plyr)
library(dplyr)
library(gtsummary)
library(multipleNCC)
library(intsurv)
library(Hmisc)
library(tidyr)
library(ggplot2)
library(rlang)
library(intsurv)
library(prodlim)
library("boot")
library(gridExtra)

# Directories: 
dir_project <- paste0("//",file.path("storage","v","vcl04","DERM","Data","ONDERZOEK","Studies","2021","EMCD21053 [StepIdent_CP_refined]","2-Data","Methodological Paper"))
dir_data <- file.path(dir_project,"Data")
dir_Lakeman_pdata <-  file.path(dir_data,"Raw Data","Lakeman_","Storage","Final database and script") # Where preprocessed Lakeman data are stored
dir_results <- file.path(dir_project,"Analyses","Results")

# dir github repository:
dir_scripts <-  file.path("C:","Users","b.rentroia","Documents","CP_refined Repository","ncc-evaluation","Code")

# Script ID
script_id <-  "RW1"

# Auxilliary functions:
source(file.path(dir_scripts,"Auxilliary Scripts","2023_03_03 Weighted performance metrics functions.R"))
source(file.path(dir_scripts,"Auxilliary Scripts","2023_03_07 NCC cohorts functions.R"))
#---------------------------------------------------

#---------------------------------------------------
# 1. Load datasets
#---------------------------------------------------
# We are going to analyze the ERGO cohort, but only a reduced subset of variables are of interest to us. 
# Lakeman et al have pre-processed this dataset, and applied the BOADICEA model to it. 
# For this real-world illustration, we load the variables of interest and the BOADICEA risk probabilities, which were merged into a single dataframe (ERGOcohort) in a previous script. 
load(file.path(dir_data,"Preprocessed Data","bc_data_with_BOADICEA_estim.RData"))

# The BOADICEA model can only be applied to subjects that are younger than 70y old and do not have breast cancer:
ERGOcohort <-  ERGOcohort[ERGOcohort$no_incident_case=="1" & ERGOcohort$startage_lower70==1,]
#---------------------------------------------------


#---------------------------------------------------
# 2. Define variables of interest
#---------------------------------------------------
outcome_num_name <-  "Affected_BC10"
outcome_FUP_name <-  "TimeToEvent_max10"
  
outcome_num <-  ERGOcohort[,outcome_num_name]
outcome_FUP <-  ERGOcohort[,outcome_FUP_name]

  
# Summarize the patient characteristics:
id_vars <-  c("ergoid","startdat","fp_date_lastcontact")
table1 <- 
  ERGOcohort [,!colnames(ERGOcohort)%in%id_vars]%>%
  mutate(Affected_BC10 = factor(Affected_BC10,levels=c("0","1"))) %>%
  sjlabelled::remove_all_labels() %>%
  mutate_if(is.factor,fct_explicit_na,na_level = "Unknown") %>% #This ensures that unknown variable is included in the tables
  tibble %>%
  tbl_summary(
    by = Affected_BC10, type = list(where(is.integer) ~ "continuous",where(is.numeric) ~ "continuous")) %>%
  add_n() %>% # add column with total number of non-missing observations
  add_p() %>% # test for a difference between groups
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels()

table1 %>%
  as_flex_table() %>%
  flextable::save_as_docx(path=file.path(dir_results,paste0(script_id," Subject characteristics stratified.docx")))
#---------------------------------------------------

#---------------------------------------------------
# 3. Calculate performance metrics of the BOADICEA model in the full cohort
#---------------------------------------------------
# Discrimination metrics:
# The following commented code gives the same results as what is reported in the publication of Lakeman 2020:  
# data.Cindex<- with(ERGOcohort, cbind(TimeToEvent_max10, Affected_BC10, BC_10_RF_PRS))
# Cindex.eval<-Est.Cval(data.Cindex,10)
# Cindex.eval$Dhat 
# Note: the function Est.Cval computes Uno's concordance metric, which is similar to Harrell's c-index but not exactly the same. 

# We compute the performance metrics for all 10-year risk estimates based on:  age only (BC_10_age), age+risk factors (BC_10_RF), PRS (BC_10_PRS), and all variables (BC_10_RF_PRS)
predictions_of_interest <- c("BC_10_age","BC_10_RF","BC_10_PRS","BC_10_RF_PRS")
performance_metrics_fullcohort <- as.data.frame(matrix(0,ncol=length(predictions_of_interest),nrow=7))
colnames(performance_metrics_fullcohort)<-predictions_of_interest
rownames(performance_metrics_fullcohort)<- c("Cindex","OEratio","CSlope","SE","SP","PPV","NPV")
performance_metrics_fullcohort_cis <-performance_metrics_fullcohort
x1<-1:nrow(ERGOcohort) # auxilliary variable to compute confidence intervals

for (pi in predictions_of_interest){
  set.seed(123) # needed for reproducibility of confidence intervals
  # C-index with bootstrap confidence interval for each metric:
  
  # C-index:
  cind_est <-round(computeCindex(outcome_FUP,outcome_num,ERGOcohort[,pi],10,weights=NULL),5)
  b1<-boot(x1,function(u,i) computeCindex(outcome_FUP[i],outcome_num[i],ERGOcohort[i,pi],10,weights=NULL),R=1000)
  b1_ci<- boot.ci(b1,type=c("norm","basic","perc"))
  performance_metrics_fullcohort_cis["Cindex",pi]<-  paste0(cind_est, " (",round(b1_ci$percent[4],5),"-",round(b1_ci$percent[5],5),")")
  performance_metrics_fullcohort["Cindex",pi]<-cind_est
  
  # Calibration metrics:
  # Observed to events ratio:
  oe_ratio_est <- round(computeOEratio(outcome_FUP,outcome_num,ERGOcohort[,pi],10),5)
  b2 <-boot(x1,function(u,i) computeOEratio(outcome_FUP[i],outcome_num[i],ERGOcohort[i,pi],10,weights=NULL),R=1000)
  b2_ci<- boot.ci(b2,type=c("norm","basic","perc"))
  performance_metrics_fullcohort_cis["OEratio",pi]<-paste0(oe_ratio_est, " (",round(b2_ci$percent[4],5),"-",round(b2_ci$percent[5],5),")")
  performance_metrics_fullcohort["OEratio",pi]<-oe_ratio_est
  
  # Calibration slope:
  cslope_est <-round(calibrationSlope(outcome_FUP,outcome_num,ERGOcohort[,pi],10),5)
  b3 <-boot(x1,function(u,i) calibrationSlope(outcome_FUP[i],outcome_num[i],ERGOcohort[i,pi],10,weights=NULL),R=1000)
  b3_ci<- boot.ci(b3,type=c("norm","basic","perc"))
  performance_metrics_fullcohort_cis["CSlope",pi]<-paste0(cslope_est, " (",round(b3_ci$percent[4],5),"-",round(b3_ci$percent[5],5),")")
  performance_metrics_fullcohort["CSlope",pi]<-cslope_est
  
  # Other discriminative metrics at cutoff 3%:
  cutoff_metrics <- round(computeDiscMetrics(outcome_FUP,outcome_num,ERGOcohort[,pi],0.03,10,weights=NULL),5)
  b4<-boot(x1,function(u,i) computeDiscMetrics(outcome_FUP[i],outcome_num[i],ERGOcohort[i,pi],0.03,10,weights=NULL),R=1000)
  b4_ci_SE<- boot.ci(b4,type=c("norm","basic","perc"),index=1)
  performance_metrics_fullcohort_cis["SE",pi]<-paste0(cutoff_metrics[1], " (",round(b4_ci_SE$percent[4],5),"-",round(b4_ci_SE$percent[5],5),")")
  
  b4_ci_SP<- boot.ci(b4,type=c("norm","basic","perc"),index=2)
  performance_metrics_fullcohort_cis["SP",pi]<-paste0(cutoff_metrics[2], " (",round(b4_ci_SP$percent[4],5),"-",round(b4_ci_SP$percent[5],5),")")
  
  b4_ci_PPV<- boot.ci(b4,type=c("norm","basic","perc"),index=3)
  performance_metrics_fullcohort_cis["PPV",pi]<-paste0(cutoff_metrics[3], " (",round(b4_ci_PPV$percent[4],5),"-",round(b4_ci_PPV$percent[5],5),")")
  
  b4_ci_NPV<- boot.ci(b4,type=c("norm","basic","perc"),index=4)
  performance_metrics_fullcohort_cis["NPV",pi]<-paste0(cutoff_metrics[4], " (",round(b4_ci_NPV$percent[4],5),"-",round(b4_ci_NPV$percent[5],5),")")
  
  performance_metrics_fullcohort[names(cutoff_metrics),pi]<-cutoff_metrics 
}
write.csv(performance_metrics_fullcohort,file.path(dir_results,paste0(script_id," Performance metrics FULL cohort.csv")))
write.csv(performance_metrics_fullcohort_cis,file.path(dir_results,paste0(script_id," Performance metrics FULL cohort CIS.csv")))

# Calibration plot on the full cohort:
pdf(file.path(dir_results,paste0(script_id," calibration plot full cohort.pdf")),width=8,height=6)
cp_fullcohort <- calib_plot(ERGOcohort[,pi],outcome_num,outcome_FUP,10,875,lim=c(0,0.15),pl=TRUE,weights=NULL,log_plot=FALSE,main_plot="",surv_or_risk = "risk",new_plot = TRUE,dots_col = "darkred",show_segments=FALSE)
dev.off()

# Decision curve analysis on the full cohort:
col_vars <-  c("BC_10_age" = "darkgoldenrod2" ,"BC_10_RF"="darkred","BC_10_RF_PRS" = "blue","BC_10_PRS" = "pink")
pdf(file.path(dir_results,paste0(script_id," dca full cohort.pdf")),width=8,height=7)
dc_fullcohort<-dca_survival(ERGOcohort,predictions_of_interest,col_vars[predictions_of_interest],outcome_num,outcome_FUP,c(-0.02,0.05),c(0,15),10,rep(1,nrow(ERGOcohort)),var_labels=c("Age","Age+RF","Age+PRS","Age+RF+PRS"))
dev.off()
#---------------------------------------------------

#---------------------------------------------------
# 4. Generate NCC cohorts from the full cohort and evaluate model:
#---------------------------------------------------
# We use a loop to derive NCC cohorts and evaluate the performance of the BOADICEA model in each of them.
# We perform some sensitivity analysis:
## number of repetitions of derivation of NCC cohorts: nrp
## number of controls per case: ncontrols
## evaluate combined model or subset of components: pr_names
## evaluate different methods for computing sampling weights: wegh_methods
# Note that we dont perform all the combinations of these parameters, to reduce the number of analyses.
# The entire loop takes around 1h to run.

# Repetitions
for (nrp in c(100,500)){
  if(nrp ==100){
    ncontrols <-c(1,2) 
  }else if (nrp==500){
    ncontrols <-1  # only test one scenario, since it takes quite some time to perform 500 repetitions.
  }
  
  for (nc in ncontrols){
    all_results <- as.data.frame(matrix(NA,nrow=0,ncol = 8))
    if(nrp ==100&nc ==1 ){
      pr_names <- c("BC_10_RF","BC_10_PRS","BC_10_RF_PRS") # For this comparison, 100 repetitions is enough.
    }else{
      pr_names <-c("BC_10_RF_PRS") 
    }
    for(p_name in pr_names){
      if(p_name %in% c("BC_10_RF","BC_10_PRS")){
          wgh_methods <-c("KM")
      }else if (p_name =="BC_10_RF_PRS"){
        if (nrp ==100 & nc ==1){
          wgh_methods <- c("KM","GAMprob","GLMprob","Prevalence","Prevalence_KM")
        }else if (nc ==2){
          wgh_methods <- "KM"
        }
      }
      
      
      # Results directory:  
      dir_results_nr_nc <- file.path(dir_results,paste0("nrep",nrp,"_nc",nc))
      if(!dir.exists(dir_results_nr_nc)){dir.create(dir_results_nr_nc)}

      for (wgh_met in wgh_methods){
        include_vars <-  c("TimeToEvent_max10",p_name,"ergoid","startage")
        
        # Each setup has its own directory: 
        dir_results_p <- file.path(dir_results_nr_nc,p_name)
        if(!dir.exists(dir_results_p)){dir.create(dir_results_p)}
        
        dir_results_wgh <- file.path(dir_results_p,wgh_met)
        if(!dir.exists(dir_results_wgh)){dir.create(dir_results_wgh)}
        
        file_name_pref_dca = file.path(dir_results_wgh,paste0(script_id," DCA"))
        file_name_pref_cp = file.path(dir_results_wgh,paste0(script_id," CPl"))
        
        # 4.1 Scenario A: Only Incidence Density Sampling
        #---------------------------------------------------
        performance_metrics_A <- evaluateModelInNCC(outcome_FUP_name,outcome_num_name,ERGOcohort,include_vars,match_vars=NULL,ncontrols=nc,pt_id="ergoid",wgh_method=wgh_met,123,predictions_name=p_name,tp=10,nrep=nrp,ylim_dca = c(-0.02,0.05),xlim_dca=c(0,15),ylim_dca_un_fact=15,file_name_dca = paste0(file_name_pref_dca," SA"),c(0,0.2),5,ylim_cp_un_fact=5,file_name_cp = paste0(file_name_pref_cp," SA"),cp_fullcohort = cp_fullcohort,dc_fullcohort=dc_fullcohort) 
        #---------------------------------------------------
        
        # 4.2 Scenario B: Incidence Density Sampling, with matching on model variable
        #---------------------------------------------------
        ERGOcohort$round_age <-round(ERGOcohort$startage, 0)
        ERGOcohort$Round_Risk_RF <- factor(round(ERGOcohort$BC_10_RF,2))
        
        #Note: many attempts were tried to match on risk factors, but correlation with the outcome was too low, therefore we decided to match on the non genetic risk estimate.
        #ERGOcohort$children_categories <-cut(ERGOcohort$Number_children, c(-Inf, 0, Inf), labels=c("None", "At least 1"))
        #ERGOcohort$OAC_or_Hormone <- factor(ifelse(ERGOcohort$OAC_ever =="Yes"|ERGOcohort$Hormone_ever =="Yes","Yes","No"),levels=c("No","Yes"))
        #ERGOcohort$Extreme_Riskfactors <- cut(ERGOcohort$BC_10_RF, c(-Inf, quantile(ERGOcohort$BC_10_RF,0.20), quantile(ERGOcohort$BC_10_RF,0.80),Inf), labels=c("<20%","20-80", ">80%"))
        
        if (p_name=="BC_10_RF"){
          # In this case we cannot match on the non genetic risk factor, so we match based on age.
          performance_metrics_B_round_RF <- evaluateModelInNCC(outcome_FUP_name,outcome_num_name,ERGOcohort,include_vars,match_vars="round_age",ncontrols=nc,pt_id="ergoid",wgh_method=wgh_met,123,predictions_name=p_name,tp=10,nrep=nrp,ylim_dca = c(-0.02,0.05),xlim_dca=c(0,15),ylim_dca_un_fact=15,file_name_dca = paste0(file_name_pref_dca," SB age"),c(0,0.2),5,ylim_cp_un_fact=5,file_name_cp = paste0(file_name_pref_cp," SB age"),cp_fullcohort = cp_fullcohort,dc_fullcohort=dc_fullcohort) 
        }else{
          performance_metrics_B_round_RF <- evaluateModelInNCC(outcome_FUP_name,outcome_num_name,ERGOcohort,include_vars,match_vars="Round_Risk_RF",ncontrols=nc,pt_id="ergoid",wgh_method=wgh_met,123,predictions_name=p_name,tp=10,nrep=nrp,ylim_dca = c(-0.02,0.05),xlim_dca=c(0,15),ylim_dca_un_fact=15,file_name_dca = paste0(file_name_pref_dca," SA RF"),c(0,0.2),5,ylim_cp_un_fact=5,file_name_cp = paste0(file_name_pref_cp," SA RF"),cp_fullcohort = cp_fullcohort,dc_fullcohort=dc_fullcohort)
        }
         
        # 4.3 Scenario C: Incidence Density Sampling, with matching on admin variable
        #---------------------------------------------------
        performance_metrics_C <- evaluateModelInNCC(outcome_FUP_name,outcome_num_name,ERGOcohort,include_vars,match_vars="Study",ncontrols=nc,pt_id="ergoid",wgh_method=wgh_met,123,predictions_name=p_name,tp=10,nrep=nrp,ylim_dca = c(-0.02,0.05),xlim_dca=c(0,15),ylim_dca_un_fact=15,file_name_dca = paste0(file_name_pref_dca," SC"),c(0,0.2),5,ylim_cp_un_fact=5,file_name_cp = paste0(file_name_pref_cp," SC"),cp_fullcohort = cp_fullcohort,dc_fullcohort=dc_fullcohort) 
        
      
        #---------------------------------------------------
        # 4.4. Summarize everything into tables/figures
        #---------------------------------------------------
        # Summarize all metrics in a single dataframe:
        results_list <-list("NCC-NM"=performance_metrics_A, "NCC-MR"=performance_metrics_B_round_RF,"NCC-MNR" = performance_metrics_C)
        results_mdf <- as.data.frame(do.call(rbind,lapply(names(results_list), function(x) cbind(pivot_longer(as.data.frame(results_list[[x]][["Performance_metrics_all"]]),      cols = everything()),x))))
        results_mdf <- results_mdf %>% mutate(Weighted = ifelse(grepl("wgh",name),"Yes","No"),
                                              Metric = revalue(gsub("wgh_|un_","",name),c("Cind"="C-index","CSlope"="Calibration slope","OE"="O/E ratio","SE"="Sensitivity","SP"="Specificity")),
                                              Full_cohort = performance_metrics_fullcohort[gsub("wgh_|un_","",name),p_name],
                                              Full_cohort_lb = as.numeric(gsub(".*\\(|-.*","",performance_metrics_fullcohort_cis[gsub("wgh_|un_","",name),p_name])),
                                              Full_cohort_ub = as.numeric(gsub(".*-|\\)","",performance_metrics_fullcohort_cis[gsub("wgh_|un_","",name),p_name])),
                                              x = factor(x,levels=c("NCC-NM","NCC-MNR","NCC-MR")),
                                              Metric = factor(Metric,levels=c("C-index","Calibration slope","O/E ratio","Sensitivity","Specificity","PPV","NPV")))%>%
                                              as.data.frame()
        # Threshold-free metrics:
        p_summary <- results_mdf %>% filter(Metric %in% c("C-index","Calibration slope","O/E ratio")) %>% ggplot(aes(x=x,y=value,col=Weighted))+geom_boxplot()+geom_hline(aes(yintercept =Full_cohort))+geom_hline(aes(yintercept =Full_cohort_lb),lty=2)+geom_hline(aes(yintercept =Full_cohort_ub),lty=2)+facet_wrap(~Metric,scales="free")+theme_bw()+xlab("Scenarios")+ylab("Metric value") 
        ggsave(file.path(dir_results_wgh,paste0(script_id," summary no thresh.pdf")),p_summary,width=10, height=3)
        
        # Threshold-based metrics:
        p_summary_wth <- results_mdf %>% filter(!Metric %in% c("C-index","Calibration slope","O/E ratio")) %>%
          ggplot(aes(x=x,y=value,col=Weighted))+geom_boxplot()+geom_hline(aes(yintercept =Full_cohort))+geom_hline(aes(yintercept =Full_cohort_lb),lty=2)+geom_hline(aes(yintercept =Full_cohort_ub),lty=2)+facet_wrap(~Metric,scales="free")+theme_bw()+xlab("Scenarios")+ylab("Metric value")
        ggsave(file.path(dir_results_wgh,paste0(script_id," summary wth thresh.pdf")),p_summary_wth,width=10, height=6)
        
        # Save results with this weighting scheme:
        results_mdf$weight_scheme <-wgh_met
        results_mdf$Model <-p_name
        all_results <- rbind(all_results,results_mdf)
        
        # Save results:
        sink(file.path(dir_results_wgh,paste0(script_id,"all metrics.txt")))
        print(results_list)
        closeAllConnections()
        #---------------------------------------------------
      }
    
    }
    
    if(nrp ==100 & nc ==1){
      
      # 5.1 Summarize all weighting schemes in one figure

      df_pl <- all_results %>%
        mutate(weight_scheme = ifelse(Weighted=="No","None",weight_scheme),
               weight_scheme = factor(weight_scheme,levels=c("None",wgh_methods))) %>%
        rename(Weights = weight_scheme)%>%
        distinct(.keep_all = TRUE) 
      
      sum_noth <-  df_pl %>% # The unweighted metrics are the same for all weighting schemes, so this line removes the repeated rows.
        filter(Metric %in% c("C-index","Calibration slope","O/E ratio"),Model =="BC_10_RF_PRS") %>%
        ggplot(aes(x=x,y=value,col=Weights))+geom_boxplot()+geom_hline(aes(yintercept =Full_cohort))+geom_hline(aes(yintercept =Full_cohort_lb),lty=2)+geom_hline(aes(yintercept =Full_cohort_ub),lty=2)+facet_wrap(~Metric,scales="free")+theme_bw()+xlab("Scenarios")+ylab("Metric value")
      
      ggsave(file.path(dir_results_nr_nc ,paste0(script_id," summary no threshold metrics all weights.pdf")),sum_noth,width=10, height=3)
      
      sum_wth <-  df_pl %>% # The unweighted metrics are the same for all weighting schemes, so this line removes the repeated rows.
        filter(!Metric %in% c("C-index","Calibration slope","O/E ratio"),Model =="BC_10_RF_PRS") %>%
        ggplot(aes(x=x,y=value,col=Weights))+geom_boxplot()+geom_hline(aes(yintercept =Full_cohort))+geom_hline(aes(yintercept =Full_cohort_lb),lty=2)+geom_hline(aes(yintercept =Full_cohort_ub),lty=2)+facet_wrap(~Metric,scales="free")+theme_bw()+xlab("Scenarios")+ylab("Metric value")
      
      ggsave(file.path(dir_results_nr_nc ,paste0(script_id," summary wth threshold metrics all weights.pdf")),sum_wth,width=10, height=6)
      
      sum_allwth <-  df_pl %>% filter(Model =="BC_10_RF_PRS") %>%# The unweighted metrics are the same for all weighting schemes, so this line removes the repeated rows.
        ggplot(aes(x=x,y=value,col=Weights))+geom_boxplot()+geom_hline(aes(yintercept =Full_cohort))+geom_hline(aes(yintercept =Full_cohort_lb),lty=2)+geom_hline(aes(yintercept =Full_cohort_ub),lty=2)+facet_wrap(~Metric,scales="free")+theme_bw()+xlab("Scenarios")+ylab("Metric value")
      
      ggsave(file.path(dir_results_nr_nc ,paste0(script_id," summary all all weights.pdf")),sum_allwth,width=10, height=6)
    
        
    
      # 5.2 Compare all prediction models:
      
      g1 <- df_pl %>% filter(Weights %in% c("None","KM"),Metric == c("C-index")) %>%# The unweighted metrics are the same for all weighting schemes, so this line removes the repeated rows.
        mutate(Model=factor(revalue(Model,c(BC_10_RF = "Age+RF",BC_10_PRS = "Age+PRS",BC_10_RF_PRS = "Age+RF+PRS")),levels =c("Age+RF","Age+PRS","Age+RF+PRS"))) %>%
        ggplot(aes(x=x,y=value,fill=Model))+geom_boxplot()+geom_hline(aes(yintercept =Full_cohort),lty=2)+facet_wrap(~Weights,ncol = 2)+theme_bw()+xlab("Scenarios")+ylab("Metric value")
      
      g2 <- df_pl %>% filter(Weights %in% c("None","KM"),Metric == c("Calibration slope")) %>%# The unweighted metrics are the same for all weighting schemes, so this line removes the repeated rows.
        mutate(Model=factor(revalue(Model,c(BC_10_RF = "Age+RF",BC_10_PRS = "Age+PRS",BC_10_RF_PRS = "Age+RF+PRS")),levels =c("Age+RF","Age+PRS","Age+RF+PRS"))) %>%
        ggplot(aes(x=x,y=value,fill=Model))+geom_boxplot()+geom_hline(aes(yintercept =Full_cohort),lty=2)+facet_wrap(~Weights,ncol = 2)+theme_bw()+xlab("Scenarios")+ylab("Metric value")
      
      g3 <- df_pl %>% filter(Weights %in% c("None","KM"),Metric == c("O/E ratio")) %>%# The unweighted metrics are the same for all weighting schemes, so this line removes the repeated rows.
        mutate(Model=factor(revalue(Model,c(BC_10_RF = "Age+RF",BC_10_PRS = "Age+PRS",BC_10_RF_PRS = "Age+RF+PRS")),levels =c("Age+RF","Age+PRS","Age+RF+PRS"))) %>%
        ggplot(aes(x=x,y=value,fill=Model))+geom_boxplot()+geom_hline(aes(yintercept =Full_cohort),lty=2)+facet_wrap(~Weights,ncol = 2)+theme_bw()+xlab("Scenarios")+ylab("Metric value")
      
      sum_all_mods <-grid.arrange(arrangeGrob(g1,top=grid::textGrob("C-index",x=0.43)), arrangeGrob(g2,top=grid::textGrob("Calibration slope",x=0.43)),arrangeGrob(g3,top=grid::textGrob("O/E ratio", x = 0.4, hjust = 0)), nrow=3)
      ggsave(file.path(dir_results_nr_nc ,paste0(script_id," comparison different models.pdf")),sum_all_mods,width=8, height=8)
    
    }
  }
}
