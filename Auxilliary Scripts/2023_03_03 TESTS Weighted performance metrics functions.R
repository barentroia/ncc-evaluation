#---------------------------------------------------
# Description: Tests of functions for the computation of performance metrics
# Author: B. Rentroia Pacheco
#---------------------------------------------------

library(survival)


# Setup:
#---------------------------------------------------
# We are going to use the ovarian dataset from the survival package:

# Calibration metrics and decision curve analyses require model predictions that can be interpreted as probabilities:
# Therefore, we generate predictions that are between 0 and 1:
ovarian$ecog.ps_01 = ovarian$ecog.ps-1
ovarian$age_01 = ovarian$age/max(ovarian$age)

# Some auxilliary variables for the decision curve analyses
vars=c("ecog.ps_01","age_01")
col_vars = c("darkslategray3","seagreen")
names(col_vars) = vars
timepoint_interest=1227

# Generate random weights for controls, case weights=1
weights_ov = rep(1,length(ovarian$fustat))
weights_ov[which(ovarian$fustat==0)]=sample(c(2,1.5,3,5,1),length(which(ovarian$fustat==0)),replace=TRUE)

# Dataset with repeated rows, based on the weights. We will use this dataset to check that weighted performance metrics are similar to the unweighted metrics in this dataset:
ovarian_rep=ovarian[rep(1:nrow(ovarian),weights_ov),]
#---------------------------------------------------



# Check Calibration metrics:
#---------------------------------------------------
# Calibration slope using log-log of survival probability, which would yield the same calibration slope as directly using the linear predictor.
# We chose this way because we can use the survival probabilities directly.
predictions = 1-ovarian$age_01-0.01 #Note we remove 0.01 to avoid infinite log-log predictions. It has no influence in the conclusions because these model probabilites are fake anyway.
lp.val <- log(-log(1-predictions))
f.val <- coxph(Surv(futime, fustat) ~ lp.val,data=ovarian)  
slope <- f.val$coefficients[1]
slope
# This is the same as with my function:
calibrationSlope(ovarian$futime,ovarian$fustat,predictions,tp=NULL,weights=NULL)

# Weighted slope is the same as providing weights to the cox model:
calibrationSlope(ovarian$futime,ovarian$fustat,predictions,tp=NULL,weights=weights_ov)
coxph(Surv(futime, fustat) ~ lp.val,data=ovarian,weights = weights_ov)$coefficients[1] 
# And it's similar to repeating observations according to weights
calibrationSlope(ovarian_rep$futime,ovarian_rep$fustat,rep(predictions,weights_ov),tp=NULL,weights=NULL)
#---------------------------------------------------

# Check Decision curve analysis:
#---------------------------------------------------

# Unweighted curves should be the same as the ones plotted with the dca function of the dcurves package:
hline_check = -0.23
hline_check_2 =0.065

p_dca = dca(Surv(futime,fustat) ~ age_01+ecog.ps_01, 
            data = ovarian,time=timepoint_interest,
            label = list(age_01="Age",ecog.ps_01="ECOG")) 

as_tibble(p_dca) %>%
  dplyr::filter(!is.na(net_benefit)) %>%
  ggplot(aes(x = threshold, y = net_benefit, color = label)) +geom_line()+
  coord_cartesian(ylim = c(-0.5,0.6)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  theme_bw() +
  theme(legend.position="none")+
  geom_hline(yintercept=hline_check ,
             color = "red", size=1)+
  geom_hline(yintercept=hline_check_2 ,
             color = "red", size=1)

dca_survival(ovarian,vars,col_vars,ovarian$fustat,ovarian$futime,c(-0.25,0.6),c(0,100),timepoint_interest,rep(1,nrow(ovarian)))
  lines(seq(0,100,1),rep(hline_check,101),col="red")
  lines(seq(0,100,1),rep(hline_check_2,101),col="red")

# The decision curve analyses are identical. The last part of the graph is different because our function extends survival probability computation when there are no observations with follow up higher than the indicated.  

# Weighted curves should be similar to the ones plotted with the dca function of the dcurves package, with repeated subjects:
# With repetition of patients based on weights

hline_check =-0.41
# Dcurves package with repeated rows:
p_dca = dca(Surv(futime,fustat) ~ age_01+ecog.ps_01, 
            data = ovarian_rep,time=timepoint_interest,
            label = list(age_01="Age",ecog.ps_01="ECOG")) 

as_tibble(p_dca) %>%
  dplyr::filter(!is.na(net_benefit)) %>%
  ggplot(aes(x = threshold, y = net_benefit, color = label)) +geom_line()+
  coord_cartesian(ylim = c(-0.5,0.6)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  theme_bw() +
  theme(legend.position="none")+
  geom_hline(yintercept=hline_check ,
             color = "red", size=1)
# DCA coded by myself:
dca_survival(ovarian_rep,vars,col_vars,ovarian_rep$fustat,ovarian_rep$futime,c(-0.5,0.6),c(0,100),timepoint_interest,rep(1,nrow(ovarian_rep)))
lines(seq(0,100,1),rep(hline_check,101),col="red")

# DCA coded by myself, using weighted approach:
dca_survival(ovarian,vars,col_vars,ovarian$fustat,ovarian$futime,c(-0.5,0.6),c(0,100),timepoint_interest,weights_ov)
lines(seq(0,100,1),rep(hline_check,101),col="red")

