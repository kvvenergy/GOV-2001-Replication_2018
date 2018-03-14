######################################################
##  Brad DeWees, Averell Schmidt, Kala Viswanathan  ##
##  Replication Paper                               ##
##  Gov 2001, Spring 2018                           ##
######################################################

library(dplyr)
library(ggplot2)
library(stargazer)
library(haven)
library(sandwich)
#library(mvtnorm)
#library(Zelig)
library(rms)
library(multiwayvcov)

#############################
############################# 
## Part 1. Replicate DGP in R
#############################
#############################

# Set Up Workspace
## Kala and Brad: You will need to change the working drive to your file path
setwd("/Users/Avery/Dropbox/GOV 2001 Replication Paper/Replication Paper/Warren_IO_2014_rep/")
## Kala's path...
## Brad's path...

samb_data <- read_dta("SambanisJCR2004_replicationdataset.dta") # Load data

############
## Correct Typos in DV Data
############

samb_data$cowcode[samb_data$cowcode ==	347] <- 345
samb_data$cowcode[samb_data$cowcode ==	364] <- 365
samb_data$cowcode[samb_data$cowcode ==	511] <- 580
samb_data$cowcode[samb_data$cowcode ==	529] <- 530
samb_data$cowcode[samb_data$cowcode ==	817] <- 816
samb_data$cowcode[samb_data$cowcode ==	818] <- 816

samb_data$atwarns_rev <- samb_data$atwarns
samb_data$atwarns_rev[samb_data$cowcode ==	710 & samb_data$year == 1957] <- 1
samb_data$atwarns_rev[samb_data$cowcode ==	710 & samb_data$year == 1958] <- 1
samb_data$atwarns_rev[samb_data$cowcode ==	710 & samb_data$year == 1959] <- 1

samb_data$warstnsb_rev <- samb_data$warstnsb
samb_data$atwarns_rev[samb_data$cowcode ==	775 & samb_data$year == 1952] <- 0

samb_data <- subset(samb_data, (is.na(samb_data$warstnsb_rev) == FALSE))

####################
## Merge with Warren's Data
####################

warren_data <- read_dta("Warren.MediaData.v1.dta") # Load data
data <- merge(warren_data, samb_data, by=c("cowcode","year"), all.x = FALSE, all.y = TRUE)

####################
## GENERATE DEPENDENT VARIABLE 
####################

data$onset <- data$warstnsb_rev

####################
## Peace Year Splines
####################

## From here I am not able to replicate Warren's code in R
## He uses a btscs package in STATA that is not the exact same in R
## We will have to do a little more work to replicate this. 
## What I did is below. 

#library(splines)
library(pltesim) # needed for btscs
data$atwar <- data$atwarns_rev
dat <- btscs(df = data, event = "atwar", t_var = "year", cs_unit = "cowcode", pad_ts = FALSE)

## Warren's code in stata does everythign in one line, I think we will need a few in R
#btscs atwar year cowcode, gen(pcyrs) nspline(3);
#rename _spline1 spline1;
#rename _spline2 spline2;
#rename _spline3 spline3;

## The following code is just me inspecting the output of the btscs command
dat1 <- tbl_df(dat)%>%  # drop unnecessary variables 
  subset(select = c(cowcode, year, spell_time)) 
dat1$t <- dat1$spell_time + 1

#############################
############################# 
## Part 2. Replicate Analysis
#############################
#############################

rm(list = ls()) # clear workspace
# import Warren's cleaned data
data <- read_dta("Warren_IO_reg_data.dta") # Load data

## Begin by replicating analysis in Warren_IO_regs.do

###############
## Reproduce Table 1 
##############

## The Catch: R doesn't have a pre-packaged function for cluster-robust standard errors in R, so I found on on the following website:
## https://stackoverflow.com/questions/33927766/logit-binomial-regression-with-clustered-standard-errors

## One possible extension is taking log of IV, this eliminates his findings
#data$mdi <- log(data$mdi)
#hist(data$mdi)
#data$mdi[data$mdi <	-100] <- NA # a handful of countries had MDIs of zero, resulting in -inf when logged

# Get coefficients using GLM function 
model1 <- glm(onset ~ lgdpl + larea + lmtn + lpopl + oil2l + deml + deml2 + ethfracl + relfracl + pcyrs + spline1 + spline2 + spline3, data = data, family=binomial(logit))
model2 <- glm(onset ~ mdi + larea + lmtn + lpopl + oil2l + deml + deml2 + ethfracl + relfracl + pcyrs + spline1 + spline2 + spline3, data = data, family = binomial(logit))
model3 <- glm(onset ~ mdi + lgdpl + larea + lmtn + lpopl + oil2l + deml + deml2 + ethfracl + relfracl + pcyrs + spline1 + spline2 + spline3, data = data, family = binomial(logit))
model4 <- glm(onset ~ mdi + lgdpl + larea + lmtn + lpopl + oil2l + deml + deml2 + ethfracl + relfracl + telephli + pcyrs + spline1 + spline2 + spline3 , data = data, family = binomial(logit))
model5 <- glm(onset ~ mdi + lgdpl + larea + lmtn + lpopl + oil2l + deml + deml2 + ethfracl + relfracl + litli + pcyrs + spline1 + spline2 + spline3, data = data, family = binomial(logit))
model6 <- glm(onset ~ mdi + lgdpl + larea + lmtn + lpopl + oil2l + deml + deml2 + ethfracl + relfracl + seduli + pcyrs + spline1 + spline2 + spline3, data = data, family = binomial(logit))
model7 <- glm(onset ~ mdi + lgdpl + larea + lmtn + lpopl + oil2l + deml + deml2 + ethfracl + relfracl +  pfl + pcyrs + spline1 + spline2 + spline3, data = data, family = binomial(logit))

## Compute cluster robust standard errors by hand using Molly Roberts code
# http://projects.iq.harvard.edu/files/gov2001/files/sesection_5.pdf
# I converted Molly's code into a function that computes SEs clustered on cowcode
library(dplyr); library(sandwich)
cluster_se <- function(dataset, mod){
  library(dplyr); library(sandwich); library(lmtest) # import relevant packages
  vars <- variable.names(mod) # getvariable from model
  vars <- subset(vars, (!vars == "(Intercept)")) # drop intercept
  vars <- c("cowcode", vars) # add cluster variable
  se_data <- tbl_df(dataset) %>% subset(select = vars) # subset data
  se_data <- na.omit(se_data) # drop rows with missing values
  countries <- as.matrix(se_data$cowcode) # extract cluster ID
  m <- length(unique(countries)) # number of clusters
  k <- length(coef(mod)) # number of IVs
  u <- estfun(mod) 
  u.clust <- matrix(NA, nrow=m, ncol=k)
  for(j in 1:k){
    u.clust[,j] <- tapply(u[,j], countries, sum)
  }
  bread <-vcov(mod) # extract variance-covariance matrix
  cl.vcov <- bread %*% ((m / (m-1)) * t(u.clust) # Make cluster robust matrix
                        %*% (u.clust)) %*% + bread
  mod <- coeftest(mod, cl.vcov) # get test coefficients
  return(mod)
}
# Apply function to models
model1 <- cluster_se(data, model1)
model2 <- cluster_se(data, model2)
model3 <- cluster_se(data, model3)
model4 <- cluster_se(data, model4)
model5 <- cluster_se(data, model5)
model6 <- cluster_se(data, model6)
model7 <- cluster_se(data, model7)

# Produce Table 1 
table1 <- stargazer(model1, model2, model3, model4, model5, model6, model7,  
                    title="Warren 2014 Table 1: Logistic Regression -- civil war onset", 
                    omit.stat = c("adj.rsq", "f"), omit = c("spline1", "spline2", "spline3"), digits = 4,
                    add.lines = list(c("Splines (1-3)", "N/S", "N/S", "N/S", "N/S", "N/S", "N/S", "N/s")),
                    notes = c("Note: All independent variables lagged by one year. Robust standard errors in parentheses.", 
                              "N/S indicates splines were included, but were not significant.", 
                              "$^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01"),
                    notes.append = FALSE, notes.align = "l")

################
## Reproduce Density Plot
################

rm(list = ls()) # clear workspace
# import Warren's cleaned data
data <- read_dta("Warren_IO_reg_data.dta") # Load data

data_1989 <- subset(data, (year == 1989)) 
# Warren's Density Plot
densplot1 <- ggplot(data_1989, aes(x = mdi)) + 
  geom_density(position="identity", alpha=0.6, stat = "density") + 
  scale_x_continuous(name = "Media Density Index, 1989 (by country)") + #, 
                     #breaks = seq(0, 340, 100),
                     #limits=c(0, .006)) +
  scale_y_continuous(name = "Kernel density") + #,
                     #breaks = seq(0, 0.006, .002),
                     #limits=c(0, .006)) +
  ggtitle("Kernel density plot") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 12))
  #guides(fill=guide_legend(title="U.S. Regions"))
densplot1
ggsave("densplot1.jpeg", plot = last_plot(), width = 6, height = 6, units = "in")
# Note this is not quite correct. Looks similar but slightly off


# Produce Kernel density for all years 
densplot2 <- ggplot(data, aes(x = mdi)) + 
  geom_density(position="identity", alpha=0.6, stat = "density") + 
  scale_x_continuous(name = "Media Density Index, all years (by country)") + #, 
  #breaks = seq(0, 340, 100),
  #limits=c(0, .006)) +
  scale_y_continuous(name = "Kernel density") + #,
  #breaks = seq(0, 0.006, .002),
  #limits=c(0, .006)) +
  ggtitle("Kernel density plot") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 12))
#guides(fill=guide_legend(title="U.S. Regions"))
densplot2
ggsave("densplot2.jpeg", plot = last_plot(), width = 6, height = 6, units = "in")


# Produce Kernel density for all years with logged MDI
data$mdi <- log(data$mdi)
data$mdi[data$mdi <	-100] <- NA # a handful of countries had MDIs of zero, resulting in -inf when logged

densplot3 <- ggplot(data, aes(x = mdi)) + 
  geom_density(position="identity", alpha=0.6, stat = "density") + 
  scale_x_continuous(name = "Logged Media Density Index, all years (by country)") + #, 
  #breaks = seq(0, 340, 100),
  #limits=c(0, .006)) +
  scale_y_continuous(name = "Kernel density") + #,
  #breaks = seq(0, 0.006, .002),
  #limits=c(0, .006)) +
  ggtitle("Kernel density plot") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 12))
#guides(fill=guide_legend(title="U.S. Regions"))
densplot3
ggsave("densplot3.jpeg", plot = last_plot(), width = 6, height = 6, units = "in")


################
## Figure 3 - ROC Curves for Models 1-3
################
#Model1
library(ROCR)
vars <- variable.names(model1)
vars <- subset(vars, (!vars == "(Intercept)")) 
data_roc1 <- data
data_roc1 <- tbl_df(data)%>%  
  subset(select = c(onset, lgdpl,larea, lmtn, lpopl, oil2l, deml, deml2, ethfracl, relfracl,pcyrs, spline1, spline2,spline3))

data_roc1 <- na.omit(data_roc1) # drop rows with missing values

predict1 <- predict(model1, type = "response")
pred1 <- prediction(predict1, data_roc1$onset)    
perf1 <- performance(pred1, measure = "tpr", x.measure = "fpr")     

#Check that AUC values match
auc_ROCR1 <- performance(pred1, measure = "auc")
auc_ROCR1 <- auc_ROCR1@y.values[[1]]
auc_ROCR1

#Model2
vars2 <- variable.names(model2)
vars2 <- subset(vars2, (!vars2 == "(Intercept)")) 
data_roc2 <- data
data_roc2 <- tbl_df(data)%>%  
  subset(select = c(onset,mdi,larea,lmtn,lpopl,oil2l,deml,deml2,ethfracl,relfracl,pcyrs,spline1, spline2,spline3))
data_roc2 <- na.omit(data_roc2) # drop rows with missing values

predict2 <- predict(model2, type = "response")
pred2 <- prediction(predict2, data_roc2$onset)    
perf2 <- performance(pred2, measure = "tpr", x.measure = "fpr")  

#Check that AUC values match
auc_ROCR2 <- performance(pred2, measure = "auc")
auc_ROCR2 <- auc_ROCR2@y.values[[1]]
auc_ROCR2

#Model3
vars3 <- variable.names(model3)
vars3 <- subset(vars3, (!vars3 == "(Intercept)")) 
data_roc3 <- data
data_roc3 <- tbl_df(data)%>%  
  subset(select = c(onset,mdi,lgdpl,larea,lmtn,lpopl,oil2l,deml,deml2,ethfracl,relfracl, spline1, spline2,spline3))

data_roc3 <- na.omit(data_roc3) # drop rows with missing values

predict3 <- predict(model3, type = "response")
pred3 <- prediction(predict3, data_roc3$onset)    
perf3 <- performance(pred3, measure = "tpr", x.measure = "fpr")     

#Check that AUC values match
auc_ROCR3 <- performance(pred3, measure = "auc")
auc_ROCR3 <- auc_ROCR3@y.values[[1]]
auc_ROCR3

#Figure 3 from the paper - plot all three ROC curves
plot(perf1, col = "grey56" )
plot(perf2, add = TRUE, col = "gray26")
plot(perf3, add = TRUE, col="gray0")
abline(0, 1) 
legend("bottomright", lwd = 2, legend=c("Controls + GDP", "Controls + Media", "Controls + GDP + Media"), col=c("grey56", "grey26","grey0"),cex=0.5,pt.cex=0.5)
ggsave("ROCplot.jpeg", plot = last_plot(), width = 6, height = 6, units = "in")

# I'm adding text here just to practice how changes get implemented in GitHub's workflow.
