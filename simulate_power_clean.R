# Supporting code for "Effectiveness of CoronaVac in the setting of high SARS-CoV-2 P.1 variant transmission in Brazil: A test-negative case-control study"
# Matt Hitchings
# No analytical power formula applies to 1:1 matched case-control study with categorical exposure
# Code is to estimate the power of a given case-control data set with specific distribution of discordant pairs through simulation

require(dplyr)
require(survival)

# Inputs: matched case control data set
# Outline: randomly assign vaccination data based on assumed odds ratios and estimate VE
assumed_ors_primary = c(1,1,0.5,0.5,0.3)
assumed_ors_atleastonedose = c(1,0.4)

# Make case-control data
run_mod = function(data,exposure,exposure_cat) {
  
  formula_unadj = formula(paste0('case~factor(',exposure,')+strata(stratum_no)'))
  
  mod_unadj = clogit(formula_unadj,data=data)
  
  ests_cis_unadj = cbind(exp(mod_unadj$coefficients),exp(confint(mod_unadj)))
  
  pvals = sapply(1:length(mod_unadj$coefficients),function(x) round(2*pnorm(-abs(mod_unadj$coefficients[x]/sqrt(vcov(mod_unadj)[x,x])),lower.tail=T),2))
  
  return(list(ests_cis_unadj,pvals))
  
}

# Read in data
# 1:1 matched case control data set, exposures are vcat_beforetest (at least one dose) and vcat_combined (1 and 2 doses)
# Variables: stratum_no, case (0,1), vcat_beforetest (0=unvaccinated, 1=at least one dose 0-13 days, 2=at least one dose 14+ days), vcat_combined (Unvaccinated, Dose 1 0-13 days, Dose 1 >=14 days, Dose 2 0-13 days, Dose 2 >=14 days)

dat$vcat_beforetest_permute = as.numeric(as.character(dat$vcat_beforetest))
dat$vcat_combined_permute = dat$vcat_combined

dat = dat %>% select(stratum_no,case,vcat_beforetest,vcat_beforetest_permute,vcat_combined,vcat_combined_permute) %>% arrange(stratum_no,case)

# Label each pair as discordant
dat$disc_atleastonedose = sapply(dat$stratum_no,
                                                function(x) 
                                                  paste(sum(dat$vcat_beforetest[dat$stratum_no==x]==0),
                                                        sum(dat$vcat_beforetest[dat$stratum_no==x]==1),
                                                        sum(dat$vcat_beforetest[dat$stratum_no==x]==2),sep='_'))

dat$disc_primary = sapply(dat$stratum_no,
                                                function(x) 
                                                  paste(sum(dat$vcat_combined[dat$stratum_no==x]=="Unvaccinated"),
                                                        sum(dat$vcat_combined[dat$stratum_no==x]=="Dose 1 0-13 days"),
                                                        sum(dat$vcat_combined[dat$stratum_no==x]=="Dose 1 >=14 days"),
                                                        sum(dat$vcat_combined[dat$stratum_no==x]=="Dose 2 0-13 days"),
                                                        sum(dat$vcat_combined[dat$stratum_no==x]=="Dose 2 >=14 days"),
                                                        sep='_'))


ves_atleastonedose = rep(NA,1000)
ps_atleastonedose = rep(NA,1000)
ves_primary = rep(NA,1000)
ves_uci_primary = rep(NA,1000)
ps_primary = rep(NA,1000)

for (i in 1:1000) {

  # 0 vs 1 pairs
  dat$vcat_beforetest_permute[dat$disc_atleastonedose=='1_1_0' & dat$case==1] = 
    rbinom(length(dat$vcat_beforetest[dat$disc_atleastonedose=='1_1_0' & dat$case==1]),1,assumed_ors_atleastonedose[1]/(1+assumed_ors_atleastonedose[1]))
  dat$vcat_beforetest_permute[dat$disc_atleastonedose=='1_1_0' & dat$case==0] = 
    sapply(dat$stratum_no[dat$disc_atleastonedose=='1_1_0' & dat$case==0],
           function(x) 1-dat$vcat_beforetest_permute[dat$stratum_no==x & dat$case==1]
    )
  
  # 0 vs 2 pairs
  dat$vcat_beforetest_permute[dat$disc_atleastonedose=='1_0_1' & dat$case==1] = 
    2*rbinom(length(dat$vcat_beforetest[dat$disc_atleastonedose=='1_0_1' & dat$case==1]),1,assumed_ors_atleastonedose[2]/(1+assumed_ors_atleastonedose[2]))
  dat$vcat_beforetest_permute[dat$disc_atleastonedose=='1_0_1' & dat$case==0] = 
    sapply(dat$stratum_no[dat$disc_atleastonedose=='1_0_1' & dat$case==0],
           function(x) 2-dat$vcat_beforetest_permute[dat$stratum_no==x & dat$case==1]
    )
  
  # 1 vs 2 pairs
  dat$vcat_beforetest_permute[dat$disc_atleastonedose=='0_1_1' & dat$case==1] = 
    1+rbinom(length(dat$vcat_beforetest[dat$disc_atleastonedose=='0_1_1' & dat$case==1]),1,assumed_ors_atleastonedose[2]/assumed_ors_atleastonedose[1]/
               (1+assumed_ors_atleastonedose[2]/assumed_ors_atleastonedose[1]))
  dat$vcat_beforetest_permute[dat$disc_atleastonedose=='0_1_1' & dat$case==0] = 
    sapply(dat$stratum_no[dat$disc_atleastonedose=='0_1_1' & dat$case==0],
           function(x) 3-dat$vcat_beforetest_permute[dat$stratum_no==x & dat$case==1]
    )
  
  t=run_mod(dat,'vcat_beforetest_permute')
  ves_atleastonedose[i] = unname(t[[1]][2,1])
  ps_atleastonedose[i] = unname(t[[2]][2])
  
  # Primary analysis
  # Unvaccinated vs Dose 1 0-13 days pairs
  vcat_labels = c("Unvaccinated","Dose 1 0-13 days","Dose 1 >=14 days","Dose 2 0-13 days","Dose 2 >=14 days")
  for (p1 in 1:4) {
    for (p2 in (p1+1):5) {
      pair = rep(0,5)
      pair[c(p1,p2)]=1
      pair = paste0(pair,collapse='_')
      dat$vcat_combined_permute[dat$disc_primary==pair & dat$case==1] = 
        sapply(1:length(dat$vcat_combined[dat$disc_primary==pair & dat$case==1]),function(x) 
          vcat_labels[c(p1,p2)][1+rbinom(1,1,assumed_ors_primary[p2]/assumed_ors_primary[p1]/(1+assumed_ors_primary[p2]/assumed_ors_primary[p1]))])
      dat$vcat_combined_permute[dat$disc_primary==pair & dat$case==0] = 
        sapply(dat$stratum_no[dat$disc_primary==pair & dat$case==0],
               function(x) setdiff(vcat_labels[c(p1,p2)],dat$vcat_combined_permute[dat$stratum_no==x & dat$case==1])
        )
    }
  }
  
  dat$vcat_combined_permute = relevel(factor(dat$vcat_combined_permute),ref='Unvaccinated')
  
  t=run_mod(dat,'vcat_combined_permute')
  ves_primary[i] =  unname(t[[1]][3,1])
  ves_uci_primary[i] =  unname(t[[1]][3,3])
  ps_primary[i] =  unname(t[[2]][3])
  
}

summary(ves_atleastonedose)
mean(ps_atleastonedose<=0.05)
summary(ves_primary)
summary(ves_uci_primary)
mean(ves_uci_primary<=0.5)
mean(ps_primary<=0.05)
