% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probanofix.R
\name{probanofix}
\alias{probanofix}
\title{Probabitility of event when the follow-up is no fixed}
\usage{
probanofix(surv,time,duraccrual,limit,gamma)
}
\arguments{
\item{surv}{Survival estimates}

\item{time}{Time estimate}

\item{duraccrual}{Accrual duration, expressed in t time units}

\item{limit}{Time limit to estimate the survival probability}

\item{gamma}{the probability of observing an event by time t}
}
\value{
probanofix: event probability at time limit
}
\description{
In the design with variable follow-up, each subject is followed until the end of the study (Tf = infinite), i.e. subjects who are 
enrolled at the beginning of the enrolment phase are followed for a longer time than subjects who are enrolled later. 
When there are no drop outs (i.e. = (0,0), the probability of failure in each arm can be directly estimated using the formulation 
proposed by K Kim and A.A. Tsiatis (1990).
}
\details{
The probanofix function estimates the probability of event when the follow-up is no fixed
}
\references{
K. Kim, A.A. Tsiatis, Study duration for clinical trials with survival response and early stopping rule, Biometrics 46 (1990) 81-92
}

