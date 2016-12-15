############################################
#               probafix.R                 #
############################################
#'
#' @title Probabitility of event when the follow-up is fixed
#'  
#' @description In a fixed follow-up design, each subject can only be followed during a fixed period (Tf < infinite) 
#' and then goes off study. Using similar reasoning to K. Kim and A.A. Tsiatis (1990), it is easy to compute the probability.
#' 
#' @usage probafix(surv,time,duraccrual,durfollow,limit,gamma)
#' 
#' @param surv Survival estimates
#' @param time Time estimate
#' @param duraccrual Accrual duration, expressed in t time units
#' @param durfollow Follow-up duration
#' @param limit Time limit to estimate the survival probability
#' @param gamma the probability of observing an event by time t
#'  
#' @details The probafix function estimates the probability of event when the follow-up is fixed. 
#' 
#' @return probafix: event probability at time limit
#' 
#' @references K. Kim, A.A. Tsiatis, Study duration for clinical trials with survival response and early stopping rule, Biometrics 46 (1990) 81-92
#' 
#' @name probafix
#' @export 
#' 
NULL

probafix <- function(surv,time,duraccrual,durfollow,limit,gamma)
{ 
	lambda <- -log(surv)/time
	limitstar <- limit-durfollow
	
	# Case 1: limit <= Accrual duration
	cond1 <- 0<=limit
	cond2 <- limit<=duraccrual
	cond3 <- limitstar<0
	cond <- cond1+cond2+cond3
	if (cond==3)
	{
		#	print("Case 1")
		part <- limit-1/(lambda+gamma)+exp(-(lambda+gamma)*limit)/(lambda+gamma)
		probafix <- part/limit
	}
	
	# Case 2: limit <= Accrual duration limit>= Follow-up duration & limit <= Accrual duration
	cond1 <- 0<=limit
	cond2 <- limit<=duraccrual
	cond3 <- 0<=limitstar
	cond4 <- limitstar<= (duraccrual-durfollow)
	cond <- cond1+cond2+cond3+cond4
	if (cond==4)
	{
		#	print("Case 2")
		part1 <- limit-1/(lambda+gamma)
		part2 <- limit-durfollow-1/(lambda+gamma)
		probafix <- exp(-(lambda+gamma)*durfollow)*part2
		probafix <- part1-probafix
		probafix <- probafix/limit
	}
	
	# Case 3: 
	cond1 <- duraccrual<limit
	cond2 <- limit<=(duraccrual+durfollow)
	cond3 <- limitstar<0
	cond <- cond1+cond2+cond3
	if (cond==3)
	{
		#	print("Case 3")
		part <- exp((lambda+gamma)*duraccrual)-1
		part <- exp(-(lambda+gamma)*limit)/(lambda+gamma)*part
		probafix <- duraccrual-part
		probafix <- probafix/duraccrual
	}
	
	# Case 4:
	cond1 <- duraccrual<limit
	cond2 <- limit<=(duraccrual+durfollow)
	cond3 <- 0<=limitstar
	cond4 <- limitstar<=duraccrual
	cond <- cond1+cond2+cond3+cond4
	if (cond==4) 
	{
		#	print("Case 4")
		part1 <- duraccrual-(limit-durfollow)*exp(-(lambda+gamma)*durfollow)
		part2 <- exp((lambda+gamma)*duraccrual)-exp((lambda+gamma)*(limit-durfollow))
		part2 <- exp(-(lambda+gamma)*limit)/(lambda+gamma)*part2
		probafix <- (part1-part2)/duraccrual
	}
	
	# Case 5
	cond1 <- (duraccrual+durfollow)<limit
	cond2 <- duraccrual<limitstar
	cond <- cond1+cond2
	if (cond==2)
	{
		print("Case 5")
		part <- 1-exp(-(lambda+gamma)*durfollow)
		probafix <- duraccrual*part/duraccrual
	}
	
	probafix <- probafix*lambda/(lambda+gamma)
	return(probafix)
}
