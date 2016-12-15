\name{coprimary-package}
\alias{coprimary-package}
\alias{coprimary}
\docType{package}
\title{Sample size calculation for two primary time-to-event endpoints in phase III clinical trials}
\description{The coprimary R package computes the required number of patients for two time-to-event endpoints as primary endpoint. 
This R package contains six functions to check the consistency of the design and computes the sample size with one or two time-to-event endpoints. 
Both endpoints can be one time to event endpoint, such as Overall Survival (OS), Progression Free Survival (PFS) or the health health-related quality of life (HRQoL), or two time-to-event 
endpoints, which could be PFS and time to HRQoL score deterioration as example.}
\details{

\tabular{ll}{
Package:\tab coprimary\cr
Type: \tab package\cr
Version:\tab 1.0\cr
Date: \tab 2016-09-20\cr
licence: \tab GPL(>=3.3.2)\cr
}
}
\author{
Alhousseiny PAM\cr 
Maintainer: <alhousseiny.pam@gmail.com>
}

\references{
Legault, Claudine. Analyzing multiple endpoints with a two-stage group sequential design in clinical trials. Diss. University of North Carolina at Chapel Hill, 1991.\cr
Chow, Shein-Chung, Hansheng Wang, and Jun Shao. Sample size calculations in clinical research. CRC press, 2007.\cr
Filleron, T., Gal, J., & Kramar, A. (2012). Designing group sequential randomized clinical trials with time to event end points using a R function. Computer methods 
and programs in biomedicine, 108(1), 113-128.
}

\examples{

#####################################################################################
############ Design superiority: two-sided with two co-primary endpoints ############
#####################################################################################
## - For endpoint 1: 2 target variables for the health related quality of life with 3-year 
## rate without HRQoL deterioration Se=0.75 and Sc=0.67, alpha1=c(0.01,0.01)
## - For endpoint 2: 4-year survival rates Se=0.86 and Sc=0.80, alpha2=c(0.015,0.015)
## - with accrual duration of 3 years, study duration of 6 years, power=0.90, 
## look=c(2,c(1,1),0.5), and default values i.e  pe=0.5, fup=0, dropout=0

nc1 <- ncoprimary(design=c(1,2),Survhyp1=c(1,5,0.75,0.67),Survhyp2=c(1,5,0.86,0.80),
alpha1=c(0.01,0.01),alpha2=c(0.015,0.015),duraccrual=3,durstudy=6,power=0.90,
look=c(2,c(1,1),0.5),dqol=2)


#####################################################################################
############ Design superiority: one-sided with two co-primary endpoints ############
#####################################################################################
## - For endpoint 1: 2-year hazard ratio hype=0.86 and Sc=0.62, alpha1=0.05
## - For endpoint 2: 3-year survival rates hype=0.81 and Sc=0.57, alpha2=0.05
## - with accrual duration of 2 years, study duration of 10 years and default values i.e
## power=0.90, pe=0.5, look=1, fup=0, dropout=0, dqol=0


nc2 <- ncoprimary(design=c(1,1),Survhyp1=c(2,2,0.86,0.62),Survhyp2=c(2,3,0.81,0.57),
alpha1=0.05,alpha2=0.05,duraccrual=2,durstudy=10)

#####################################################################################
###############   Design non-inferiority with one primary endpoint ##################
#####################################################################################
## 5-year rate without HRQoL deterioration are equal under the alternative hypothesis, 
## i.e Se=0.60 and Sc=SeA=0.70, with alpha=0.05, accrual duration of 4 years,study duration 
## of 8 years, two interim analysis after the occurence 1/3 and 2/3 of the events, 3 target 
## variables for the health related quality of life and default values i.e power=0.80, pe=0.5, 
## fup=0, dropout=0

ns <- nsurvival(design=c(2),Survhyp=c(1,5,0.60,0.70, 0.70),alpha=0.05,duraccrual=4,
durstudy=8,look=c(3,c(1,1),c(1/3,2/3)), dqol=3)

}

\keyword{sample size}
\keyword{multiple endpoints} 
\keyword{co-primary}
\keyword{clinical trial}
