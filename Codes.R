library(adjustedCurves)
library(boot)
library(broom)
library("geepack")
library(here)
library(MatchIt)
library(tableone)
library(sjPlot)
library(survey)

# IMPORT DATASET
                                        
rhc <- as.data.frame(read.csv("rhc.csv", header=T))
rhc$A <- ifelse(rhc$swang1 =="No RHC", 0, 1)

print("Variables with missing data - ")
names(which(colSums(is.na(rhc))>0))

print("Proportion of Missing data - ")
p <- ((colMeans(is.na(rhc)))*100); p[which(p> 0)]

rhc$surv <- (rhc$dthdte - rhc$sadmdte)+1
#If date of death missing, then the patient is alive by the end of the follow-up.
rhc[c("surv")][is.na(rhc[c("surv")])] <- 9999 

rhc$surv30 <- ifelse(rhc$surv <32, 0, 1)
rhc$surv60 <- ifelse(rhc$surv <61, 0, 1)
rhc$surv180 <- ifelse(rhc$surv <180, 0, 1)
rhc$d30 <- 1- rhc$surv30; rhc$d60 <- 1- rhc$surv60; rhc$d180 <- 1- rhc$surv180 
                                      
covs <- c("surv30","surv60","surv180","death")
factors <- c("surv30","surv60","surv180")
table <- CreateTableOne(vars=covs, factorVars = factors, strata="swang1", data=rhc, test=TRUE)
print("Crude outcome - ")
table



                                        
##PROPENSITY SCORE - CODE
d <- glm(A ~ age +factor(sex) +factor(race) +edu +factor(income) +factor(ninsclas) +factor(cat1) +resp +card 
         +neuro +gastr +renal +meta +hema +seps +trauma +ortho +das2d3pc +factor(dnr1) +factor(ca) +surv2md1 +aps1 
         +scoma1 +wtkilo1 +temp1 +meanbp1 + resp1 +hrt1 + pafi1 +paco21 +ph1 +wblc1 +hema1 +sod1 + pot1+ crea1 
         +bili1 +alb1 +factor(cardiohx) +factor(chfhx) +factor(dementhx) +factor(psychhx) +factor(chrpulhx) 
         +factor(renalhx) +factor(liverhx) +factor(gibledhx) +factor(malighx) +factor(immunhx) +factor(transhx) 
         +factor(amihx), data=rhc, family=binomial(link=logit))

deno <- predict(d, type="response")
hist(deno[rhc$A==1], xlim=c(0,1), ylim=c(0, 225), col="red", breaks=40, main='Propensity scores', xlab="")
hist(deno[rhc$A==0], add=T, col=rgb(0, 0.1, 1, 0.5), breaks= 40)
legend('topright', legend=c('RHC','No RHC'), lwd=2, col=c('red','blue'))

                                        


## IPTW - code
  #numerator 
n <- glm(A ~ 1, data=rhc, family=binomial(link=logit)); nume <- predict(n, type="response")
                                      
rhc$SIPTW <- (rhc$A==1) * nume/deno + (rhc$A==0) * (1-nume)/(1-deno)
summary(rhc$SIPTW)
                                      
  # Weighted covariate balance:
covs <- c("age", "sex", "race", "edu", "income", "ninsclas", "cat1", "resp", "card", 
          "neuro", "gastr", "renal", "meta", "hema", "seps", "trauma", "ortho", 
          "das2d3pc", "dnr1", "ca", "surv2md1", "aps1", "scoma1", "wtkilo1", "temp1", 
          "meanbp1", "resp1", "hrt1", "pafi1", "paco21", "ph1", "wblc1", "hema1", 
          "sod1",  "pot1", "crea1", "bili1", "alb1", "cardiohx", "chfhx", "dementhx", 
          "psychhx", "chrpulhx", "renalhx", "liverhx", "gibledhx", "malighx", "immunhx", "transhx", "amihx")

weighted <- svydesign(ids=~0, data=rhc, weights=rhc$SIPTW)
table.w <- svyCreateTableOne(vars=covs, strata="A", data=weighted, smd=TRUE, test=F)
print(table.w, smd=TRUE)

 #Compare to Table 1 of the *JAMA* paper.
                                      
  ### IPTW - Weighting the outcomes
MSM <- geeglm(surv30 ~ swang1, family=binomial("log"), data=rhc,
              weights=SIPTW, std.err = 'san.se', id=ptid, corstr="independence")
#NOTE: the log-binomial regression doesn't always converge. 
#In practice we use a logistic regression for a binary outcomes.
#Then we extract predicted probabilities to estimates the RR, or RD.
#Look at lines  96 to 120

MSM.1 <- geeglm(surv30 ~ swang1, family=binomial("logit"), data=rhc,
                weights=SIPTW, std.err = 'san.se', id=ptid, corstr="independence")
tab_model(MSM, MSM.1) 
                                      
MSM.60 <- geeglm(surv60 ~ swang1, family=binomial("log"), data=rhc,
                 weights=SIPTW, std.err = 'san.se', id=ptid, corstr="independence")
MSM.180 <- geeglm(surv180 ~ swang1, family=binomial("log"), data=rhc,
                  weights=SIPTW, std.err = 'san.se', id=ptid, corstr="independence")
tab_model(MSM, MSM.60, MSM.180) 
                                      
#CREATE A COPY OF THE DATASET WHERE EVERYONE IS EXPOSED
A1  <- rhc; A1$swang1 <- "RHC"  
#CREATE A COPY OF THE DATASET WHERE EVERYONE IS NON-EXPOSED
A0  <- rhc; A0$swang1 <- "No RHC"
#EXTRACT THE OUTCOME PROBABILITIES
P.A1 <- predict(MSM.1, newdata=data.frame(A1), type="response")
P.A0 <- predict(MSM.1, newdata=data.frame(A0), type="response") 
#COMPUTE THE MEASURE OF EFFECT OF INTEREST
RD30 <- mean(P.A1) - mean(P.A0)
RR.30 <- mean(P.A1) / mean(P.A0);
OR.30 <- mean(P.A1/(1-P.A1)) / mean(P.A0/(1-P.A0))
OR.30a <- mean(P.A0/(1-P.A0)) / mean(P.A1/(1-P.A1))
round(cbind(RD30, RR.30, OR.30, OR.30a), 2)
                                      
 
#BOOTSTRAP IS USED FOR CONFIDENCE INTERVALS                                     
MSM.b = function(data,indices)
  {dat=data[indices,]
  Num.Mod <- glm(A ~ 1, data=dat, family=binomial(link=logit))
  dat$num <- predict(Num.Mod, type="response")
  PSA.Mod <- glm(A ~ age +factor(sex) +factor(race) +edu +factor(income) 
                 +factor(ninsclas) +factor(cat1) +resp +card +neuro +gastr +renal 
                 +meta +hema +seps +trauma +ortho +das2d3pc +factor(dnr1) 
                 +factor(ca) +surv2md1 + aps1 +scoma1 +wtkilo1 +temp1  +meanbp1 
                 +resp1 +hrt1 + pafi1 +paco21 +ph1 +wblc1 +hema1 +sod1 + pot1+ crea1
                 +bili1 +alb1 +factor(cardiohx) +factor(chfhx) +factor(dementhx) 
                 +factor(psychhx) +factor(chrpulhx) +factor(renalhx) +factor(liverhx) 
                 +factor(gibledhx) +factor(malighx) +factor(immunhx) +factor(transhx) 
                 +factor(amihx), data=dat, family=binomial(link=logit))
  dat$den     <- predict(PSA.Mod, type="response")
  dat$IPTW <- ifelse(dat$A==1, dat$num/dat$den, ifelse(dat$A==0, (1-dat$num)/(1-dat$den), NA))
  MSM <- geeglm(surv30 ~ A, family=binomial, data=dat,
                weights=IPTW, std.err = 'san.se', id=ptid, corstr="independence")
  MSM.ATE = exp(as.numeric(coefficients(MSM)[2]))
  }
   
# Get original estimate, by plugging in indices 1:n
 MSM.b(rhc, indices =1:nrow(rhc))
# Draw 100 bootstrap sample estimates
 boot.out = boot(rhc, MSM.b, 100)
# compute confidence intervals using percentile method
 boot.ci(boot.out, type = "perc", conf = 0.95)

 
 
 
## Matching
hist(deno[rhc$A==1], xlim=c(0,1), ylim=c(0, 225), col="red", breaks=40, main='Propensity scores', xlab="")
hist(deno[rhc$A==0], add=T, col=rgb(0, 0.1, 1, 0.5), breaks= 40)
legend('topright', legend=c('RHC=1','RHC=0'), lwd=2, col=c('red','blue'))

library(MatchIt)
formula <- (A ~ age +factor(sex) +factor(race) +edu +factor(income) +factor(ninsclas) +resp +card 
            +neuro +gastr +renal +meta +hema +seps +trauma +ortho +das2d3pc +factor(dnr1) +factor(ca) 
            +surv2md1 +aps1 +scoma1 +wtkilo1 +temp1 +meanbp1 + resp1 +hrt1 + pafi1 +paco21 +ph1 +wblc1 
            +hema1 +sod1 +pot1 +crea1 +bili1 +alb1 +factor(cardiohx) +factor(chfhx) +factor(dementhx) 
            +factor(psychhx) +factor(chrpulhx) +factor(renalhx) +factor(liverhx) +factor(gibledhx) 
            +factor(malighx) +factor(immunhx) +factor(transhx) +factor(amihx))
m.out1 <- matchit(formula, data = rhc, method = "nearest", exact= ~ cat1, distance = "glm", caliper = .15)
                                      
#Check the Instructions of the MatchIt paxkages here: [MatchIt](https://cran.r-project.org/web/packages/MatchIt/vignettes/MatchIt.html)
                                      

# Checking balance after NN matching
plot(m.out1, type = "jitter", interactive = FALSE)
                                      
plot(m.out1, type = "density", interactive = FALSE,
     which.xs = ~age + hrt1) #For some covariates

plot(summary(m.out1))
m.out1

Match <- match.data(m.out1)
PSM.30 <- geeglm(surv30 ~ swang1, family=binomial("log"), data=Match,
                 weights=weights, std.err = 'san.se', id=subclass, corstr="independence")
PSM.60 <- geeglm(surv60 ~ swang1, family=binomial("log"), data=Match,
                 weights=weights, std.err = 'san.se', id=subclass, corstr="independence")
PSM.180 <- geeglm(surv180 ~ swang1, family=binomial("log"), data=Match,
                  weights=weights, std.err = 'san.se', id=subclass, corstr="independence")
tab_model(PSM.30, PSM.60, PSM.180) 
                                      

# DIRECT STANDARDISATION
  #1. Fit the outcome model
OutM <- glm(surv30 ~ swang1 + age +factor(sex) +factor(race) +edu +factor(income) 
            +factor(ninsclas) +factor(cat1) +resp +card +neuro +gastr +renal +meta +hema 
            +seps +trauma +ortho +das2d3pc +factor(dnr1) +factor(ca) +surv2md1 +aps1 +scoma1
            +wtkilo1 +temp1 +meanbp1 + resp1 +hrt1 + pafi1 +paco21 +ph1 +wblc1 +hema1 +sod1 
            + pot1 +crea1 +bili1 +alb1 +factor(cardiohx) +factor(chfhx) +factor(dementhx) 
            +factor(psychhx) +factor(chrpulhx) +factor(renalhx) +factor(liverhx) 
            +factor(gibledhx) +factor(malighx) +factor(immunhx) +factor(transhx) 
            +factor(amihx), data=rhc, family=binomial(link=logit))


  #2. Averaging the exposure effect over the covariate distribution of the standard population**
Risk_A1 <- predict(OutM, newdata=A1, type="response") 
Risk_A0 <- predict(OutM, newdata=A0, type="response") 
RD.30ds <- mean(Risk_A1) - mean(Risk_A0)
RR.30ds <- mean(Risk_A1) / mean(Risk_A0)

round(cbind(RD30, RD.30ds, RR.30, RR.30ds), 3)

#SURVIVAL ANALYSES  
library(adjustedCurves)
rhc$group <- as.factor(rhc$swang1); rhc$event <- ifelse(rhc$death=="Yes",1,0)  
surv <- adjustedsurv(data=rhc,
                        variable="group",
                        ev_time="surv",
                        event="event",
                        method="iptw_km",
                        treatment_model=rhc$SIPTW,
                        conf_int=TRUE)

plot(surv, conf_int=TRUE, linetype=TRUE, max_t=181)

# library(psycModel)
# html_to_pdf(file_path = "Slides.html")
# html_to_pdf(dir = "Applied Example")

