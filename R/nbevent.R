
############################################
#               nbevent.R                  #
############################################
#'
#' @title Number of events estimates 
#' 
#' @description To determine the sample size N in clinical trials with time to event endpoint, it is necessary to proceed in two steps. In the first step, 
#' the numbers of events that need to be observed (e) are computed. In the second step, we determine the number of patients necessary to observe the number of events required.
#' This function computes the number of event for one-time-to event.
#' 
#' @usage nbevent(hypsurv,pe,alfa,beta,design)
#' 
#' @param hypsurv For Superiority=c(Sc,Se); for Non inferiority=c(Sc,Se,SeA); for Equivalence=c(Sc, Se),
#'  with Sc is survival rate in the control arm; Se is survival rate in experimental arm; SeA is the survival rate in the experimental arm under the alternative hypothesis. 
#' @param pe Proportion (ratio) of patients assigned to the experimental arm (with 0<pe<1). 
#' @param alfa Type I error, for Non inferiority, Equivalence and 1-sided superiority, alfa is a vector of length one. 
#' For 2-sided superiority, alfa is a vector to length two c(alpha.low, alpha.up).
#' @param beta Probability of a type II error.
#' @param design Superiority=c(1,sided)[with sided=1 if 1-sided and 2 if 2-sided]; Non inferiority=c(2); Equivalence=c(1,1) 
#' 
#' @details The nbevent function computes the required number of events to determine the number of patients.
#' 
#' @return E: Number of events
#' @return h: Hazard Ratio under null hypothesis(HR=log(Se)/log(Sc))
#' @return h.alt: Hazard Ratio under alternative hypothesis (h.alt=log(SeA)/log(Sc))
#' 
#' 
#' @name nbevent
#' @export 
#' @importFrom stats qnorm

#' @references Chow, S. C., Shao, J., Wang, H. (2003). Sample Size Calculation in Clinical Research. New York: Marcel Dekker.
#' @references Schoenfeld. Sample-size formula for the proportional-hazards regression model. Biometrics. 1983 39<499>503.


nbevent <- function(hypsurv,pe,alfa,beta,design)
{

	Sc <- hypsurv[1]
	Se <- hypsurv[2]
	if (design[1]==2)
	{SeA <- hypsurv[3]}
	
	
	# Trial type / one-sided or two-sided
	if (design[1]==2)
	{sided <- 1}
	else 
	{sided <- design[2]}
	
	# Quantile Alfa and Beta	/ Hazar Ratio
	palfa <- qnorm(1-alfa/sided)
	pbeta <- qnorm(1-beta)
	h <- log(Se)/log(Sc)
	if(design[1]==2)
	{h.alt <- log(SeA)/log(Sc)}
	
	# Number of events
	## Superiority
	if (design[1]==1)
	{
		part1 <- palfa+pbeta
		part2 <- log(h)*(pe*(1-pe))^(1/2)
		E <- (part1*1/part2)^2
	}
	
	## Non Inferiority
	else 
	{	
		part1 <- palfa+pbeta
		part2 <- pe*(1-pe)*(log(h)-log(h.alt))^2
		E <- part1^2/part2
	}
	if (design[1]==1)
	{return(list(E=E,h=h))}
	else 
	{return(list(E=E,h=h,h.alt=h.alt))}
}


