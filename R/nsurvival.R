############################################
#               nsurvival.R                #
############################################
#'                      
#' @title Sample size calculation in clinical trials with one primary survival endpoint 
#' @description 
#' nsurvival() is used to determine the sample size for one time to event endpoint, such as Overall Survival (OS), Progression Free Survival
#' or the health related quality of life (HRQoL). If it is HRQoL, several HRQoL dimension can be considered.
#' @usage 
#' nsurvival(design,Survhyp,alpha,duraccrual,durstudy,power,pe,look,fup,dropout,dqol)
#' 
#' 
#' @param design Superiority=c(1,sided)[with sided=1 if 1-sided and 2 if 2-sided]; Non inferiority=c(2); Equivalence=c(3) 
#' @param Survhyp For Superiority=c(thyp,t,hype,Sc); for Non inferiority=c(thyp,t,hype,Sc,hypeA); for Equivalence=c(t,delta,Sc): parameters at time t 
#' if thyp=1 then hype is survival rate in experimental arm under the null hypothesis and hypeA is the survival rate in the experimental arm under the alternative hypothesis; 
#' if thyp=2 then hype is the hazard ratio under the null hypothesis and hypeA is the hazard ratio under the alternative hypothesis; Sc is survival rate in the control arm; 
#' delta is the log hazard ratio equivalence margin. When endpoint is HRQoL, the survival rate is replaced by the rate of patients without HRQoL deterioration.
#' @param alpha Type I error, for Non inferiority, Equivalence and 1-sided superiority is a vector of length one. 
#' For 2-sided superiority is a vector to length two c(alpha.low, alpha.up). 
#' @param duraccrual Accrual duration, expressed in number of days, months or years
#' @param durstudy Study duration, expressed in number of days, months or years
#' @param power 1- Probability of a type II error. Default value=0.80.
#' @param pe Proportion (ratio) of patients assigned to the experimental arm (with 0<pe<1). Default value = 0.5.
#' @param look The number of interim analyses, c(1) for one final analysis; c(nb, bound, timing) for at least one interim analyses with bound=c(bound.eff,bound.fut):1-sided or bound=c(bound.lown,bound.up):2-sided.
#' nb the number of planned looks, bound.eff and bound.fut corresponds to the type of boundaries used for efficacy (i.e. reject H0) and futility (i.e. reject H1). 
#' bound.fut=0: No futility monitoring, 1: Lan deMets O.Brien Fleming, 2: Lan deMets Pocock. 
#' bound.low and bound.up the type of lower and upper boundaries used (1: Lan deMets O.Brien Fleming, 2: Lan deMets Pocock, 3: O.Brien Fleming, 4: Pocock). Default value = 1.
#' @param fup Follow-up information, No fixed:c(0) (follow-up until the end of study); 
#' Fixed:c(1, durfollow) with durfollow is the duration of follow-up. Default value = 0.
#' @param dropout Drop out information, No drop out=c(0); Drop out=c(1,gammae,gammac) with gammae the hazard drop out rates in experimental arm 
#' and control arm respectively. Default value = 0.
#' @param dqol number of targeted dimensions for the health related quality of life. Default value = 0.
#' 
#' @details The nsurvival function computes the sample size for one time to event endpoint, such as OS, PFS or HRQoL. 
#' HRQoL has become increasingly important in clinical trials over the past two decades.
#' 
#' @return Event: number of events estimated
#' @return Total: number of patients
#' @return Ne: number for experimental arm for each endpoint 
#' @return Nc: number for control arm for each endpoint 
#' @return HR: Hazard ratio for each endpoint
#' 
#' @examples 
#' #############################################################
#' ###############  Design superiority:one-sided ###############
#' #############################################################
#' ## 7-year survival rates Se=0.57 and Sc=0.53, alpha=0.05, accrual duration of 4 years, 
#' ## study duration of 8 years and default values i.e power=0.80, pe=0.5, look=1, fup=0, 
#' ## dropout=0, dqol=0
#'  
#' ns1 <- nsurvival(design=c(1,1),Survhyp=c(1,7,0.57,0.53),alpha=0.05,duraccrual=4,durstudy=8)
#' 
#' ############################################################
#' ############### Design superiority:two-sided ###############
#' ############################################################
#' ## 5-year rate without HRQoL deterioration Se=0.75 and Sc=0.65, alpha=c(0.04,0.01), accrual 
#' ## duration of 2 years, study duration of 6 years, power=0.90, pe=0.55, follow-up 5 years, 
#' ## 3 target variables for health related quality of life and default values i.e look=1, dropout=0
#' 
#' ns2 <- nsurvival(design=c(1,2),Survhyp=c(1,5,0.75,0.65),alpha=c(0.04,0.01),duraccrual=2,
#' durstudy=6,power=0.90,pe=0.55,fup=c(1,5),dqol=3)  
#' 
#' ###########################################################
#' ###############   Design non-inferiority ##################
#' ###########################################################
#' ## 5-year survival rates are equal under the alternative hypothesis, i.e Se=0.60 and Sc=SeA=0.70, 
#' ## with alpha=0.05, accrual duration of 4 years, study duration of 8 years, two interim analysis 
#' ## after the occurence 1/3 and 2/3 of the events and default values i.e power=0.80, pe=0.5, fup=0, 
#' ## dropout=0, dqol=0 
#' 
#' ns3 <- nsurvival(design=c(2),Survhyp=c(1,5,0.60,0.70, 0.70),alpha=0.05,duraccrual=4,
#' durstudy=8,look=c(3,c(1,1),c(1/3,2/3)))
#' 
#' ##########################################################
#' ###############    Design superiority   ##################
#' ##########################################################
#' ## 3-year rate without HRQoL deterioration Sc=0.80 and log hazard equivalence margin delta=0.1 
#' ## with alpha=0.10, accrual duration of 3 years, study duration of 5 years, drop out hazard rate 
#' ## of 0.05 per arm, 2 target variables for health related quality of life and default values i.e 
#' ## power=0.80, pe=0.5, look=1, fup=0
#' 
#' ns4 <- nsurvival(design=c(3),Survhyp=c(3,0.10,0.80),alpha=0.10,duraccrual=3,durstudy=5,
#' dropout=c(1,0.05,0.05),dqol=2)
#' 
#' @name nsurvival
#' @export 
#' 
#' @import gsDesign graphics digest plyr proto
#' @importFrom stats pnorm 


# nsurvival <- function(design,Survhyp,pe,alpha,power,duraccrual,durstudy,look,fup,dropout,dqol=NULL)
nsurvival <- function(design,Survhyp,alpha,duraccrual,durstudy,power=0.80,pe=0.5,look=1,fup=0,dropout=0,dqol=0)
{		
	beta=1-power
	
	if (dqol!=0) 
	{
		if (design[1]==1 & design[2]==2)
		{
			alfa<- c(round(alpha[1]/dqol,3),round(alpha[2]/dqol,3))
		}
		else
		{
			alfa<- round(alpha/dqol,3)			
		}
	}
	else 
	{
		alfa=round(alpha,3)
	}
	
	
	# Survival decomposition
	if (design[1]==3) {Survhyp <- c(2,Survhyp)}
	time <- Survhyp[2]
	Sc <- Survhyp[4]
	if (Survhyp[1]==1)
	{
		Se <- Survhyp[3]
		if (design[1]!=1)
		{SeA <-  Survhyp[5]}
	}
	if (Survhyp[1]==2)
	{
		time <- Survhyp[2]
		if (design[1]==1 | design[1]==2) {Se <- exp(Survhyp[3]*log(Sc))}
		if (design[1]==3) {Se <- exp(exp(-Survhyp[3])*log(Sc))}
		if (design[1]==2)
		{SeA <- exp(Survhyp[5]*log(Sc))}
	}
	if (design[1]==1)
	{Survhyp.temp <- c(Survhyp[1],time,Se,Sc)}
	if (design[1]==2)
	{Survhyp.temp <- c(Survhyp[1],time,Se,Sc,SeA)}
	if (design[1]==3)
	{Survhyp.temp <- c(time,Se,Sc)}
	
	
	if (dropout[1]==0)
	{dropout <-rep(0,3)}
	
	# Number of event / Hazard Ratio
	if (design[1]==1)
	{
		if (design[2]==1)
		{temp <- nbevent(c(Sc,Se),pe,alfa,beta,design)}
		if (design[2]==2)
		{temp <- nbevent(c(Sc,Se),pe,alfa[2],beta,c(1,1))}
		nbevent <- temp$E
		h <- temp$h
	}
	
	# Non inferiority: HR Alternative Hypothesis
	if (design[1]==2)
	{
		temp.alt <- nbevent(c(Sc,Se,SeA),pe,alfa,beta,design)
		h.alt <- temp.alt$h.alt
		nbevent <- temp.alt$E
		h <- temp.alt$h		
	}
	
	# Equivalence: HR Alternative Hypothesis
	if (design[1]==3)
	{
		temp <- nbevent(c(Sc,Se),pe,alfa,beta/2,c(1,1))
		nbevent <- temp$E
		h <- temp$h		
	}
	
	if (look[1]!=1)
	{
		if (design[1]==2)
		{design <- c(design,1)}
		
		interim.bound <- matrix(0,nrow=look[1],ncol=5+design[2])
		interim.bound[,1] <- c(look[-c(1,2,3)],1)
		
		# Si non inferiority or Superiority One Sided
		## Only H0
		if (design[2]==1 & look[3]==0)
		{
			if (look[2]==1)
			{sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa,beta=beta,sfu=sfLDOF,test.type=1)}
			if (look[2]==2)
			{sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa,beta=beta,sfu=sfLDPocock,test.type=1)}
			if (look[2]==3)
			{sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa,beta=beta,sfu="OF",test.type=1)}
			if (look[2]==4)
			{sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa,beta=beta,sfu="Pocock",test.type=1)}
			
			nbevent <- nbevent*max(sequential$n.I)
			interim.bound[,2] <- c(look[-c(1,2,3)],1)*nbevent
			interim.bound[,3] <- pnorm(-sequential$upper$bound)*design[2]
			interim.bound[,4] <- sequential$upper$bound
			if (design[1]==2 & h<1)
			{interim.bound[,4] <- -sequential$upper$bound}
			dimnames(interim.bound) <- list(c(1:look[1]),c("Information","Events","pvalue reject H0","Boundary reject H0","Analysis Time Under H0","Analysis Time Under H1"))
		}
		
		if (design[2]==1 & look[3]!=0)
		{
			if (look[2]==1)
			{
				if (look[3]==1)
				{sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa,beta=beta,sfu=sfLDOF,sfl=sfLDOF,test.type=4)}
				if (look[3]==2)
				{sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa,beta=beta,sfu=sfLDOF,sfl=sfLDPocock,test.type=4)}
			}
			
			if (look[2]==2)
			{
				if (look[3]==1)
				{sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa,beta=beta,sfu=sfLDPocock,sfl=sfLDOF,test.type=4)}
				if (look[3]==2)
				{sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa,beta=beta,sfu=sfLDPocock,sfl=sfLDPocock,test.type=4)}
			}
			nbevent <- nbevent*max(sequential$n.I)
			interim.bound[,2] <- c(look[-c(1,2,3)],1)*nbevent
			interim.bound[,3] <- sequential$upper$bound
			interim.bound[,4] <- sequential$lower$bound
			if (design[1]==2 & h<1)
			{
				interim.bound[,3] <- -sequential$upper$bound
				interim.bound[,4] <- -sequential$lower$bound
			}
			dimnames(interim.bound) <- list(c(1:look[1]),c("Information","Events","Boundary Reject H0","Boundary reject H1","Analysis Time Under H0","Analysis Time Under H1"))
		}
		
		
		if (design[2]==2)
		{
			# Inflation Parameter
			if (look[3]==1)
			{
				if (look[2]==1)
				{
					sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],astar=alfa[1],beta=beta,sfu=sfLDOF,sfl=sfLDOF,test.type=5)
					nbevent <- max(gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],beta=beta,sfu=sfLDOF,test.type=1)$n.I)*nbevent
				}
				if (look[2]==2)
				{
					sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],astar=alfa[1],beta=beta,sfu=sfLDOF,sfl=sfLDPocock,test.type=5)
					nbevent <- max(gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],beta=beta,sfu=sfLDOF,test.type=1)$n.I)*nbevent
				}
				if (look[2]==3)
				{
					sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],beta=beta,sfu=sfLDOF,sfl="OF",test.type=2)
					nbevent <- max(sequential$n.I)*nbevent
				}
				if (look[2]==4)
				{
					sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],beta=beta,sfu=sfLDOF,sfl="Pocock",test.type=2)
					nbevent <- max(sequential$n.I)*nbevent
				}
			}
			
			if (look[3]==2)
			{
				if (look[2]==1)
				{
					sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],astar=alfa[1],beta=beta,sfu=sfLDPocock,sfl=sfLDOF,test.type=5)
					nbevent <- max(gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],beta=beta,sfu=sfLDPocock,test.type=1)$n.I)*nbevent
				}
				if (look[2]==2)
				{
					sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],astar=alfa[1],beta=beta,sfu=sfLDPocock,sfl=sfLDPocock,test.type=5)
					nbevent <- max(gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],beta=beta,sfu=sfLDPocock,test.type=1)$n.I)*nbevent
				}
				if (look[2]==3)
				{
					sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],beta=beta,sfu=sfLDPocock,sfl="OF",test.type=2)
					nbevent <- max(sequential$n.I)*nbevent
				}
				if (look[2]==4)
				{
					sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],beta=beta,sfu=sfLDPocock,sfl="Pocock",test.type=2)
					nbevent <- max(sequential$n.I)*nbevent
				}
				
			}
			
			if (look[3]==3)
			{
				if (look[2]==1)
				{
					sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],beta=beta,sfu="OF",sfl=sfLDOF,test.type=2)
					nbevent <- max(sequential$n.I)*nbevent
				}
				if (look[2]==2)
				{
					sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],beta=beta,sfu="OF",sfl=sfLDPocock,test.type=2)
					nbevent <- max(sequential$n.I)*nbevent
				}
				if (look[2]==3)
				{
					sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],beta=beta,sfu="OF",sfl="OF",test.type=2)
					nbevent <- max(sequential$n.I)*nbevent
				}
				if (look[2]==4)
				{
					sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],beta=beta,sfu="OF",sfl="Pocock",test.type=2)
					nbevent <- max(sequential$n.I)*nbevent
				}
			}
			
			if (look[3]==4)
			{
				if (look[2]==1)
				{
					sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],beta=beta,sfu="Pocock",sfl=sfLDOF,test.type=2)
					nbevent <- max(sequential$n.I)*nbevent
				}
				if (look[2]==2)
				{
					sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],beta=beta,sfu="Pocock",sfl=sfLDPocock,test.type=2)
					nbevent <- max(sequential$n.I)*nbevent
				}
				if (look[2]==3)
				{
					sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],beta=beta,sfu="Pocock",sfl="OF",test.type=2)
					nbevent <- max(sequential$n.I)*nbevent
				}
				if (look[2]==4)
				{
					sequential <- gsDesign(k=look[1],timing=look[-c(1,2,3)],alpha=alfa[2],beta=beta,sfu="Pocock",sfl="Pocock",test.type=2)
					nbevent <- max(sequential$n.I)*nbevent
				}
			}
			
			interim.bound[,2] <- c(look[-c(1,2,3)],1)*nbevent
			interim.bound[,4] <-  sequential$lower$bound
			interim.bound[,5] <-  sequential$upper$bound
			dimnames(interim.bound) <- list(c(1:look[1]),c("Information","Events","pvalue reject H0","Boundary reject H0-","Boundary reject H0+","Analysis Time Under H0","Analysis Time Under H1"))
			if (alfa[1]==alfa[2])
			{interim.bound[,3] <- pnorm(-sequential$upper$bound)*2}
		}
	}
	
	# Number of patients included
	## No Fixed Follow-up
	if (fup[1]==0)
	{
		if (design[1]==1)
		{
			probae <- probanofix(Se,time,duraccrual,durstudy,dropout[2])
			probac <- probanofix(Sc,time,duraccrual,durstudy,dropout[3])
		}
		if (design[1]==2)
		{
			probae <- probanofix(SeA,time,duraccrual,durstudy,dropout[2])
			probac <- probanofix(Sc,time,duraccrual,durstudy,dropout[3])
		}
		if (design[1]==3)
		{
			probae <- probanofix(Se,time,duraccrual,durstudy,dropout[3])
			probac <- probanofix(Sc,time,duraccrual,durstudy,dropout[3])
		}
	}
	else
	## Fixed Follow-Up
	{
		if (design[1]==1)
		{
			probae <- probafix(Se,time,duraccrual,fup[2],durstudy,dropout[2])
			probac <- probafix(Sc,time,duraccrual,fup[2],durstudy,dropout[3])
		}
		if (design[1]==2)
		{
			probae <- probafix(SeA,time,duraccrual,fup[2],durstudy,dropout[2])
			probac <- probafix(Sc,time,duraccrual,fup[2],durstudy,dropout[3])
		}
		if (design[1]==3)
		{
			probae <- probafix(Se,time,duraccrual,fup[2],durstudy,dropout[2])
			probac <- probafix(Sc,time,duraccrual,fup[2],durstudy,dropout[3])
		}
		
	}
	proba <- pe*probae+(1-pe)*probac
	ntot <- nbevent/proba
	ne.temp <- ceiling(ntot*pe)
	nc.temp <- ceiling(ntot*(1-pe))
	ntot <- ne.temp+nc.temp
	
	if (fup[1]==0)
	{delaigraph <- seq(0,durstudy,by=.01)}
	else
	{delaigraph <- seq(0,max(durstudy,duraccrual+fup[2]),by=.01)}
	
	
	
	
	#############################################################################################################
	#############################################################################################################
	###################                                                              ############################
	###################        list of parameters returned                           ############################
	###################                                                              ############################
	#############################################################################################################
	#############################################################################################################

	if(fup[1]!=0)
	{
		if(dropout[1]!=0)
		{	
			if (look[1]!=1)
			{
				if (design[1]==1) 
				{
					if (design[2]==2)
					{
					l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup[2],durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
							Blower=alfa[1], Bupper=alfa[2],time=time, Sc=Sc,Se=round(Se,3), HR=round(h,3),events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe),
							subjectC=round(ntot*(1-pe)), interim=round(interim.bound,3))
					}
					if (design[2]==1)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup[2],durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,Se=round(Se,3), HR=round(h,3),events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe))
								, interim=round(interim.bound,3))
					}
				}
				
				else
				{						
					
					if (design[1]==2)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup[2],durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,Se=round(Se,3), SeA=round(SeA,3), HR=round(h,2), alt.HR=round(h.alt,3),events=round(nbevent),subjects=round(ntot), 
								subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)), interim=round(interim.bound,3))
					}
					if (design[1]==3)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup[2],durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,logHR=Survhyp[3],H01=round(exp(Survhyp[3])),H11=round(exp(Survhyp[3]),3), H02=round(exp(-Survhyp[3]),3), H12=round(exp(-Survhyp[3]),3), 
								events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)), interim=round(interim.bound,3))
					}

				}
			}
			else
			{	
				if (design[1]==1) 
				{
					if (design[2]==2)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup[2],durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								Blower=alfa[1], Bupper=alfa[2],time=time, Sc=Sc,Se=round(Se,3), HR=round(h,3),events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe),
								subjectC=round(ntot*(1-pe)))
					}
					if (design[2]==1)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup[2],durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,Se=round(Se,3), HR=round(h,3),events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)))
					}
				}
				
				else
				{						
					
					if (design[1]==2)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup[2],durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,Se=round(Se,3), SeA=round(SeA,3), HR=round(h,2), alt.HR=round(h.alt,3),events=round(nbevent),subjects=round(ntot), 
								subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)))
					}
					if (design[1]==3)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup[2],durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,logHR=Survhyp[3],H01=round(exp(Survhyp[3])),H11=round(exp(Survhyp[3]),3), H02=round(exp(-Survhyp[3]),3), H12=round(exp(-Survhyp[3]),3), 
								events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)))
					}
				}
			}
			
		}
		else
		{	
			if (look[1]!=1)
			{	
				if (design[1]==1)
				{
					if (design[2]==2)
					{
					l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup[2],durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
							Blower=alfa[1], Bupper=alfa[2],time=time, Sc=Sc,Se=round(Se,3), HR=round(h,3),events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe),
							subjectC=round(ntot*(1-pe)), interim=round(interim.bound,3))
					}
					if (design[2]==1)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup[2],durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,Se=round(Se,3), HR=round(h,3),events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe))
								, interim=round(interim.bound,3))
					}
					
				}
				if (design[1]==2)
				{
					l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup[2],durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
							time=time, Sc=Sc,Se=round(Se,3), SeA=round(SeA,3), HR=round(h,2), alt.HR=round(h.alt,3),events=round(nbevent),subjects=round(ntot), 
							subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)), interim=round(interim.bound,3))
				}
				if (design[1]==3)
				{
					l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup[2],durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
							time=time, Sc=Sc,logHR=Survhyp[3],H01=round(exp(Survhyp[3])),H11=round(exp(Survhyp[3]),3), H02=round(exp(-Survhyp[3]),3), H12=round(exp(-Survhyp[3]),3), 
							events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)), interim=round(interim.bound,3))
				}
			}
			else
			{	
				if (design[1]==1) 
				{
					if (design[2]==2)
					{
					l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup[2],durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
							Blower=alfa[1], Bupper=alfa[2],time=time, Sc=Sc,Se=round(Se,3), HR=round(h,3),events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe),
							subjectC=round(ntot*(1-pe)))
					}
					if (design[2]==1)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup[2],durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,Se=round(Se,3), HR=round(h,3),events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)))
					}
				}
				
				else
				{						
					if (design[1]==2)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup[2],durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,Se=round(Se,3), SeA=round(SeA,3), HR=round(h,2), alt.HR=round(h.alt,3),events=round(nbevent),subjects=round(ntot), 
								subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)))
					}
					if (design[1]==3)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup[2],durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,logHR=Survhyp[3],H01=round(exp(Survhyp[3])),H11=round(exp(Survhyp[3]),3), H02=round(exp(-Survhyp[3]),3), H12=round(exp(-Survhyp[3]),3), 
								events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)))
					}

				}
			}
		}
	}
	
	else
	{
		if(dropout[1]!=0)
		{
			
			if (look[1]!=1)
			{				
				if (design[1]==1) 
				{
					if (design[2]==2)
					{
					l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup,durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
							Blower=alfa[1], Bupper=alfa[2],time=time, Sc=Sc,Se=round(Se,3), HR=round(h,3),events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe),
							subjectC=round(ntot*(1-pe)), interim=round(interim.bound,3))
					}
					if (design[2]==1)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup,durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,Se=round(Se,3), HR=round(h,3),events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe))
								, interim=round(interim.bound,3))
					}
				}
				
				else
				{						
					if (design[1]==2)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup,durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,Se=round(Se,3), SeA=round(SeA,3), HR=round(h,2), alt.HR=round(h.alt,3),events=round(nbevent),subjects=round(ntot), 
								subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)), interim=round(interim.bound,3))
					}
					if (design[1]==3)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup,durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,logHR=Survhyp[3],H01=round(exp(Survhyp[3])),H11=round(exp(Survhyp[3]),3), H02=round(exp(-Survhyp[3]),3), H12=round(exp(-Survhyp[3]),3), 
								events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)), interim=round(interim.bound,3))
					}

				}
			}
			else
			{				
				if (design[1]==1) 
				{
					if (design[2]==2)
					{
					l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup,durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
							Blower=alfa[1], Bupper=alfa[2],time=time, Sc=Sc,Se=round(Se,3), HR=round(h,3),events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe),
							subjectC=round(ntot*(1-pe)))
					}
					if (design[2]==1)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup,durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,Se=round(Se,3), HR=round(h,3),events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)))
					}
				}
				
				else
				{						
					if (design[1]==2)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup,durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,Se=round(Se,3), SeA=round(SeA,3), HR=round(h,2), alt.HR=round(h.alt,3),events=round(nbevent),subjects=round(ntot), 
								subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)))
					}
					if (design[1]==3)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup,durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,logHR=Survhyp[3],H01=round(exp(Survhyp[3])),H11=round(exp(Survhyp[3]),3), H02=round(exp(-Survhyp[3]),3), H12=round(exp(-Survhyp[3]),3), 
								events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)))
					}
				}
			}
		}
		else
		{				
			if (look[1]!=1)
			{
				if (design[1]==1) 
				{
					if (design[2]==2)
					{
					l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup,durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
							Blower=alfa[1], Bupper=alfa[2],time=time, Sc=Sc,Se=round(Se,3), HR=round(h,3),events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe),
							subjectC=round(ntot*(1-pe)), interim=round(interim.bound,3))
					}
					if (design[2]==1)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup,durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,Se=round(Se,3), HR=round(h,3),events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe))
								, interim=round(interim.bound,3))
					}
				}
				
				else
				{						
					if (design[1]==2)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup,durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,Se=round(Se,3), SeA=round(SeA,3), HR=round(h,2), alt.HR=round(h.alt,3),events=round(nbevent),subjects=round(ntot), 
								subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)), interim=round(interim.bound,3))
					}
					if (design[1]==3)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup,durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,logHR=Survhyp[3],H01=round(exp(Survhyp[3])),H11=round(exp(Survhyp[3]),3), H02=round(exp(-Survhyp[3]),3), H12=round(exp(-Survhyp[3]),3), 
								events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)), interim=round(interim.bound,3))
					}
				}
			}
			else
			{				
				if (design[1]==1) 
				{
					if (design[2]==2)
					{
					l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup,durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
							Blower=alfa[1], Bupper=alfa[2],time=time, Sc=Sc,Se=round(Se,3), HR=round(h,3),events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe),
							subjectC=round(ntot*(1-pe)))
					}
					if (design[2]==1)
					{
						
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup,durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,Se=round(Se,3), HR=round(h,3),events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)))
					}
				}
				
				else
				{					
					if (design[1]==2)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup,durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,Se=round(Se,3), SeA=round(SeA,3), HR=round(h,2), alt.HR=round(h.alt,3),events=round(nbevent),subjects=round(ntot), 
								subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)))
					}
					if (design[1]==3)
					{
						l9 <- list(alpha=alfa,power=1-beta, accrual=duraccrual,followup=fup,durstuty=durstudy,ExpArm=pe,dropout=c(dropout[2],dropout[3]),look=look[1],
								time=time, Sc=Sc,logHR=Survhyp[3],H01=round(exp(Survhyp[3])),H11=round(exp(Survhyp[3]),3), H02=round(exp(-Survhyp[3]),3), H12=round(exp(-Survhyp[3]),3), 
								events=round(nbevent),subjects=round(ntot), subjectE=round(ntot*pe), subjectC=round(ntot*(1-pe)))
					}
				}
			}
		}
		
	}
	
	
	return(l9)
	
}
	

