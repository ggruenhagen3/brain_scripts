library(ggeffects)
library(sjmisc)
library(lme4)
library(splines)
# program to calculate statistical power in single nuclei studies

#set the power values to zero at the start
power_logistic<-0
power_mixed_logistic<-0
power_fisher<-0
power_ttest<-0


#PARAMATERS

#number of simulations
sim_num<-10

#total number of cells examined
cells<-10000

#mean proportion of cells displaying the transcript in group 1
p1<-.05
#standard deviation of the proportion of cells displaying the transcript among individuals in group 1
std_p1<-.02

#number of individuals per group
indiv<-10

#fold change in proportion of cells displaying the transcript between groups 1 and 2
dimorph<-2

#calculates mean proportion of cells displaying the transcript in group 2
p2<-p1*dimorph


# SIMULATION

# for loop for running the simulation
for (j in 1:sim_num)
{
  
  # making the matrix "sim" to hold the simulated values for a single run of the simulation, rows will be the number of cells
  #first column will be the animal ID 1-10, second column the outcome of whether the cell was positive or not for the transcript in group 1
  #third comlun will be the outcome of whether the cell was positive or not for the transcript in group 2 
  #this section also adds sum of total numbers of positive cells per individual  into a separate matrix named "out"
  
  sim<-matrix(NA,cells,3)
  
  #making the matrix to hold the summaries per individual used for the t-test
  out<-matrix(NA,indiv,3)
  out[,1]<-seq(1,indiv, by=1)
  
  #making the simulated values and putting them in the matrix, along with naming individuals 1-10 per group
  for (i in 1:indiv)
  {
    
    sim[ (( (cells/indiv)*(i-1) ) +1):((cells/indiv)*i),1] <-rep(1,cells/indiv)*i
    sim[ (( (cells/indiv)*(i-1) ) +1):((cells/indiv)*i),2] <-rbinom(cells/indiv,1,rnorm(1,p1,std_p1))
    sim[ (( (cells/indiv)*(i-1) ) +1):((cells/indiv)*i),3] <-rbinom(cells/indiv,1,rnorm(1,p2,std_p1))
    
    out[i,2]<-sum( sim[ (( (cells/indiv)*(i-1) ) +1):((cells/indiv)*i),2]) 
    out[i,3]<-sum( sim[ (( (cells/indiv)*(i-1) ) +1):((cells/indiv)*i),3]) 
    
    
  }
  
  
  #this section of code simply stacks the sim matrix into one that has individual ids 1-20, so individuals are in rows in stead of columns
  # it also makes the matrix into a data frame so it can be manipulated easier for running the logistic regression
  
  out3<-matrix(NA,cells*2,3)
  
  out3[(1:cells),1]<-sim[,1]
  out3[(cells+1):(cells*2),1]<-sim[,1]+indiv
  
  out3[(1:cells),2]<-sim[,2]
  out3[(cells+1):(cells*2),2]<-sim[,3]
  
  out3[(1:cells),3]<-rep(1,cells)
  out3[(cells+1):(cells*2),3]<-rep(2,cells)
  
  colnames(out3)<- c("ID", "cell", "group")
  out3<-as.data.frame(out3)
  out3$ID<-as.factor(out3$ID)
  out3$group<-as.factor(out3$group)
  
  #normal logistic regression model
  fm1<-glm(cell~group, out3, family=binomial)
  
  #keeps track of whether the logistic regression displayed a significant p-value or not
  logistic<-summary(fm1)$coefficients[2,4]
  if (logistic<0.05) {power_logistic=power_logistic+1}
  
  #mixed effects  logistic regression model entering individual as a random effect
  
  fm2 <- glmer(
    cell ~ group + (1 | ID), 
    data = out3, 
    family = binomial(link = "logit")
  )
  
  #keeps track of whether the logistic regression displayed a significant p-value or not
  
  mixed_logistic<-summary(fm2)$coefficients[2,4]
  if (mixed_logistic<0.05) {power_mixed_logistic=power_mixed_logistic+1}
  
  
  #creates a contingency table for fisher exact test analysis
  x<-matrix(NA,2,2)
  x[1,1]<-sum(sim[,2])
  x[1,2]<-cells-sum(sim[,2])
  x[2,1]<-sum(sim[,3])
  x[2,2]<-cells-sum(sim[,3])
  
  
  #fisher exact test 
  fisher<-fisher.test(x)$p.value
  
  #keeps track of whether the fisher exact test displayed a significant p-value or not
  
  if (fisher<0.05) {power_fisher=power_fisher+1}
  
  #t-test on proportions per individual
  ttest<-t.test(out[,2],out[,3])$p.value
  
  #keeps track of whether the fisher exact test displayed a significant p-value or not
  
  if (ttest<0.05) {power_ttest=power_ttest+1}
  
}


#calculates the statistical power based on the simulations
power_logistic<-power_logistic/sim_num
power_mixed_logistic<-power_mixed_logistic/sim_num
power_fisher<-power_fisher/sim_num
power_ttest<-power_ttest/sim_num

#outputs statistical power for each analysis
power_logistic
power_mixed_logistic
power_fisher
power_ttest


#-------------------------------------Not used----------------------------------------------------
#out2<-matrix(NA,indiv*2,3)

#out2[,1]<-seq(1,(indiv*2), by=1)

#out2[1:indiv,2]<-rep(1,indiv)
#out2[1:indiv,3]<-out[,2]

#out2[(indiv+1):(indiv*2),2] <- rep(2,indiv)
#out2[(indiv+1):(indiv*2),3]<-out[,3]

