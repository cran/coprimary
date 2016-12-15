############################################
#               probanofix.R               #
############################################
#'
#' @title Probabitility of event when the follow-up is no fixed
#' 
#' @description In the design with variable follow-up, each subject is followed until the end of the study (Tf = infinite), i.e. subjects who are 
#' enrolled at the beginning of the enrolment phase are followed for a longer time than subjects who are enrolled later. 
#' When there are no drop outs (i.e. = (0,0), the probability of failure in each arm can be directly estimated using the formulation 
#' proposed by K Kim and A.A. Tsiatis (1990). 
#' 
#' @usage probanofix(surv,time,duraccrual,limit,gamma)
#' 
#' @param surv Survival estimates
#' @param time Time estimate
#' @param duraccrual Accrual duration, expressed in t time units
#' @param limit Time limit to estimate the survival probability
#' @param gamma the probability of observing an event by time t
#' 
#' @details The probanofix function estimates the probability of event when the follow-up is no fixed
#' 
#' @return probanofix: event probability at time limit
#' 
#' @references K. Kim, A.A. Tsiatis, Study duration for clinical trials with survival response and early stopping rule, Biometrics 46 (1990) 81-92
#' 
#' @name probanofix
#' @export 
#' 


probanofix <- function(surv,time,duraccrual,limit,gamma)
{
	lambda <- -log(surv)/time
	
	# If time after the end of accrual
	if (limit>duraccrual)
	{
		part1 <- exp(-(lambda+gamma)*limit)/(lambda+gamma)
		part2 <- exp((lambda+gamma)*duraccrual)-1
		probanofix <- 1-part1*part2/duraccrual
	}
	else
	# If time before the end of accrual
	{
		part <- limit-1/(lambda+gamma)+exp(-(lambda+gamma)*limit)/(lambda+gamma)
		probanofix <- part/limit
	}
	
	probanofix <- probanofix*lambda/(lambda+gamma)
	return(probanofix)
}