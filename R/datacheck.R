############################################
#               datacheck.R                #
############################################

#' 
#' @title check the consistency of the parameters 
#' 
#' @description this function check the parameters required to calcul the sample size for nsurvival and ncoprimary functions. 
#' datacheck is a simple utility for carrying out parameter checks and reporting on problems or errors. 
#' 
#' @usage datacheck(design,Survhyp,pe,alfa,beta,duraccrual,durstudy,look,followup,dropout)
#' 
#' @param design Superiority=c(1,sided)[with sided=1 if 1-sided and 2 if 2-sided]; Non inferiority=c(2); Equivalence=c(3)
#' @param Survhyp For Superiority=c(thyp,t,hype,Sc); for Non inferiority=c(thyp,t,hype,Sc,hypeA); for Equivalence=c(t,delta,Sc): parameters at time t 
#' if thyp=1 then hype is survival rate in experimental arm under the null hypothesis and hypeA is the survival rate in the experimental arm under the alternative hypothesis; 
#' if thyp=2 then hype is the hazard ratio under the null hypothesis and hypeA is the hazard ratio under the alternative hypothesis; Sc is the survival rate in the control arm; 
#' delta is the log hazard ratio equivalence margin. When endpoint is HRQoL, the survival rate is replaced by the rate of patients without the HRQoL deterioration.
#' @param pe Proportion (ratio) of patients assigned to the experimental arm (with 0<pe<1)
#' @param alfa Type I error, for Non inferiority, Equivalence and 1-sided superiority, alfa is a vector of length one. 
#' For 2-sided superiority, alfa is a vector to length two c(alpha.low, alpha.up).
#' @param beta Probability of a type II error. 
#' @param duraccrual Accrual duration, expressed in number of days, months or years
#' @param durstudy Study duration, expressed in number of days, months or years
#' @param look The number interim analyses, c(1) for one final analysis; c(nb, bound, timing) for at least one interim analyses with bound=c(bound.eff,bound.fut):1-sided or bound=c(bound.low,bound.up):2-sided.
#' nb the number of planned looks, bound.eff and bound.fut corresponds to the type of boundaries used for efficacy (i.e. reject H0) and futility (i.e. reject H1). 
#' bound.fut=0: No futility monitoring, 1: Lan deMets O.Brien Fleming, 2: Lan deMets Pocock. 
#' bound.low and bound.up the type of lower and upper boundaries used (1: Lan deMets O.Brien Fleming, 2: Lan deMets Pocock, 3: O.Brien Fleming, 4: Pocock). Default value = 1.
#' @param followup Follow-up information, No fixed:c(0) (follow-up until the end of study); 
#' Fixed:c(1, durfollow) with durfollow is the duration of follow-up
#' @param dropout Drop out information, No drop out:c(0); Drop out:c(1,gammae,gammac) with gammae the hazard drop out rates in experimental arm and control arm respectively.
#' 
#' @details the datacheck function performs consistency checks on the arguments
#' 
#' @name datacheck
#' @export 


	datacheck <- function(design,Survhyp,pe,alfa,beta,duraccrual,durstudy,look,followup,dropout)
	{

		followup
		# check to all parameters
		## Design
		### Non Inferiority and superiority
		if (design[1]!=1 & design[1]!=2  & design[1]!=3)
		{stop("Design: c(1,sided) superiority, c(2) non inferiority, c(3) Equivalence ")}
		
		if (design[1]!=2 & design[1]!=3 & length(design)==1)
		{stop("Design: c(1,sided) superiority [sided=1 one-sided / sided=2 two-sided], c(2) non inferiority, c(3) equivalence")}
		
		if (design[1]!=1 & length(design)!=1)
		{stop("Design: c(1,sided) superiority [sided=1 one-sided / sided=2 two-sided], c(2) non inferiority, c(3) equivalence")}
		
		if (design[1]==1 & sum(design[2]==c(1,2))==0)
		{stop("Design: Superiority c(1,sided) with [sided=1 one-sided / sided=2 two-sided]")}
		
		## Survhyp:
		### superioty is a vector to length four
		if (design[1]==1 & length(Survhyp)!=4)
		{stop("Superiority Trial: Survival Hypothesis need to a vector of length 4: c(thyp,time,Se,Sc)")}
		
		### Non inferiority is a vector to length five
		if (design[1]==2 & length(Survhyp)!=5)
		{stop("Non Inferiority Trial: Survival Hypothesis need to a vector of length 5: c(thyp,time,Se,Sc,SeA)")}
		
		### Equivalence is a vector to length three
		if (design[1]==3) {Survhyp <- c(2,Survhyp)}
		if (design[1]==3 & length(Survhyp)!=4)
		{stop("Equivalence Trial: Survival Hypothesis need to a vector of length 3: c(time,delta,Sc)")}
		
		## hypothesis testing
		if (Survhyp[1]!=1 & Survhyp[1]!=2)
		{stop("Type of Survival Hypothesis (thyp): Survhyp=c(thyp,...) need to be: 1=Survival Rate, 2=Hazard Ratio")}
		
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
		
		## Hypothesis, Control arm (0<Sc<1)
		if (Sc<=0 | Sc>=1)
		{stop("Survival probability (Sc) must be number between 0 and 1")}
		
		## Hypothesis, Experimental arm (0<Se<1)
		if (Se<=0 | Se>=1)
		{stop("Survival probability (Se) must be number between 0 and 1")}
		
		## If non inferiority of experimental arm (SeA)
		if (design[1]==2)
		{
			if (SeA<=0 | SeA>=1)
			{stop("Survival probability (SeA) must be number between 0 and 1")}
		}
		
		## time estimates
		if (time<=0)
		{stop("time of Estimation (time) need to be positive")}
		
		## Superiority but Se<=Sc
		if (design[1]==1  & Se<=Sc)
		{stop("Superiority Trial: Surv. Experimental arm (Se) need to be > Surv. Control Arm (Sc)")}
		
		if (design[1]==2)
		{
			h.alt <- round(log(SeA)/log(Sc),2)
			h.null <- round(log(Se)/log(Sc),2)
			temp <- paste("(HRnull=",h.null,", HR Alt=",h.alt,")",sep="")
			if (h.null>1 & h.alt>=h.null)
			{
				stop(paste("Non-Inferiority Trial: if HR (Null)>1 then HR (Alt.)<HR (Null)",temp,sep=""))
			}
			
			if (h.null<1 & h.alt<=h.null)
			{
				stop(paste("Non-Inferiority Trial: if HR (Null)<1 then HR (Alt.)>HR (Null)",temp,sep=""))
			}
		}
		
		## Proportion of patient included in experimental arm
		if (pe<=0 | pe>=1)
		{stop("Proportion of patients included in the Experimental arm (pe) must be number between 0 and 1")}		
		
		## Alfa
		## Non inferiority
		if (length(alfa)!=1 & design[1]==2)
		{stop("Non inferiority Trial: Error of type I (alfa) must be number between 0 and 1")}
		
		## Equivalence
		if (length(alfa)!=1 & design[1]==3)
		{stop("Equivalence Trial: Error of type I (alfa) must be number between 0 and 1")}
		
		## Superiority One Sided
		if (length(alfa)!=1 & design[1]==1 & design[2]==1)
		{stop("Superiority Trial One-sided: Error of type I (alfa) must be number between 0 and 1")}
		
		## Alfa between 0 and 1
		if (sum(alfa<=0)!=0 | sum(alfa>=1)!=0)
		{stop("Error of type I (alfa) must be number between 0 and 1")}
		
		## Alfa between 0 and 1
		if (sum(sum(alfa)<=0)!=0 | sum(sum(alfa)>=1)!=0)
		{stop("Error of type I (alfa) must be number between 0 and 1")}
		
		## Alfa between 0 and 1
		if (length(alfa)!=2 & design[1]==1 & design[2]==2)
		{stop("Superiority Trial Two-sided: Error of type I alfa=c(alfa.low, alfa.upp) alfa.low and alfa.up must be number between 0 and 1")}
		
		## Alfa between 0 and 1
		if (alfa[1]!=alfa[2] & design[1]==1 & design[2]==2 & look[1]==1)
		{stop("Superiority Trial Two-sided without interim analysis: Error of type I alfa=c(alfa.low, alfa.upp) alfa.low need to be equal to alfa.up")}
		
		## Beta
		if (beta<=0 | beta>=1)
		{stop("Error of type II (beta) must be number between 0 and 1")}
		
		## Nb of planned Look
		if (look[1]!=1 & design[1]==3)
		{stop("Number of Planned Look : For equivalence trial, no interim analysis will be planned only a Final Analysis c(1)")}
		
		### if nb of planned look <=0
		if (look[1]<=0 & design[1]!=3)
		{stop("Number of Planned Analysis : Final Analysis c(1), Interim Analysis c(Nb Analysis,bound,timing) with Nb Analysis>0")}
		
		## length of vector
		if (length(look)!=look[1]+2 & design[1]!=3 & look[1]!=1)
		{stop("Interim Analysis: look vector of length Nb Analysis +2 c(Nb Analysis,bound,timing)")}
		
		## if superiority Two-Sided
		if (design[1]==1 & design[2]==2 & look[1]!=1 & sum(look[3]==c(1,2,3,4))!=1 & alfa[1]==alfa[2])
		{stop("Boundaries Interim Analysis : Interim Analysis c(Nb Look,bound,timing,tinterim) with bound=c(bound.low,bound.up) with bound.low and bound.up (1=Lan deMets OF/ 2= Lan deMets Pocock / 3=OF / 4=Pocock)")}
		
		if (design[1]==1 & design[2]==2 & look[1]!=1 & sum(look[3]==c(1,2))!=1 & alfa[1]!=alfa[2])
		{stop("Asymmetric Boundaries Interim Analysis : Interim Analysis c(Nb Look,bound,timing,tinterim) with bound=c(bound.low,bound.up) with bound.low and bound.up (1=Lan deMets OF/ 2= Lan deMets Pocock)")}
		
		## if superiority Two-Sided
		if (design[1]==1 & design[2]==2 & look[1]!=1 & sum(look[2]==c(1,2,3,4))!=1 & alfa[1]==alfa[2])
		{stop("Boundaries Interim Analysis : Interim Analysis c(Nb Look,bound,timing,tinterim) with bound=c(bound.low,bound.up) with bound.low and bound.up (1=Lan deMets OF/ 2= Lan deMets Pocock / 3=OF / 4=Pocock)")}
		
		if (design[1]==1 & design[2]==2 & look[1]!=1 & sum(look[2]==c(1,2))!=1 & alfa[1]!=alfa[2])
		{stop("Asymmetric Boundaries Interim Analysis : Interim Analysis c(Nb Look,bound,timing,tinterim) with bound=c(bound.low,bound.up) with bound.low and bound.up (1=Lan deMets OF/ 2= Lan deMets Pocock)")}
		
		## if superiority One Sided
		### terminal efficiency
		if (design[1]==1 & design[2]==1 & look[1]!=1 & sum(look[2]==c(1,2,3,4))!=1 & look[3]==0)
		{stop("Efficacy Boundary Interim Analysis : Interim Analysis c(Nb Look,bound,timing,tinterim) with bound=c(bound.eff,bound.fut) with bound.eff (1=Lan deMets OF/ 2= Lan deMets Pocock / 3=OF / 4=Pocock)")}
		
		if (design[1]==1 & design[2]==1 & look[1]!=1 & sum(look[2]==c(1,2))!=1 & look[3]!=0)
		{stop("Efficacy Boundary Interim Analysis : Interim Analysis c(Nb Look,bound,timing,tinterim) with bound=c(bound.eff,bound.fut) with bound.eff (1=Lan deMets OF/ 2= Lan deMets Pocock)")}
		
		### terminal futility
		if (design[1]==1 & design[2]==1 & look[1]!=1 & sum(look[3]==c(0,1,2))!=1)
		{stop("Futility Boundary Interim Analysis : Interim Analysis c(Nb Look,bound,timing,tinterim) with bound=c(bound.eff,bound.fut) with bound.fut (0= No Futility / 1=Lan deMets OF/ 2= Lan deMets Pocock)")}
		
		## if non inferiority One Sided
		### terminal efficiency
		if (design[1]==2 & look[1]!=1 & sum(look[2]==c(1,2,3,4))!=1 & look[3]==0)
		{stop("Efficacy Boundary Interim Analysis : Interim Analysis c(Nb Look,bound,timing,tinterim) with bound=c(bound.eff,bound.fut) with bound.eff (1=Lan deMets OF/ 2= Lan deMets Pocock / 3=OF / 4=Pocock)")}
		
		if (design[1]==2 & look[1]!=1 & sum(look[2]==c(1,2))!=1 & look[3]!=0)
		{stop("Efficacy Boundary Interim Analysis : Interim Analysis c(Nb Look,bound,timing,tinterim) with bound=c(bound.eff,bound.fut) with bound.eff (1=Lan deMets OF/ 2= Lan deMets Pocock)")}
		
		### terminal futility
		if (design[1]==2 & look[1]!=1 & sum(look[3]==c(0,1,2))!=1)
		{stop("Futility Boundary Interim Analysis : Interim Analysis c(Nb Look,bound,timing,tinterim) with bound=c(bound.eff,bound.fut) with bound.fut (0= No Futility / 1=Lan deMets OF/ 2= Lan deMets Pocock)")}
		
		## Si Nb Look>1 & length==1
		if (look[1]!=1 & sum(look[-c(1,2,3)]<=0)!=0)
		{stop("Timing of Interim Anbalysis : c(Nb Look,bound,timing) with 0<timing[1]<=..<=timing[nblook-1]<1")}
		
		## if Nb Look>1 & length==1
		if (look[1]!=1 & sum(look[-c(1,2,3)]>=1)!=0)
		{stop("Timing of Interim Anbalysis : c(Nb Look,bound,timing) with 0<timing[1]<=..<=timing[nblook-1]<1")}
		
		## Follow-up
		### Follow-up[1]==0 | 1
		if (sum(followup[1]==c(0,1))==0)
		{stop("Parameter followup: c(0) No Fixed Follow-up, c(1,durfup) Fixed follow with duration durfup (durfup>0)")}
		
		### No Fixed Follow-up 
		if (followup[1]==0 & length(followup)!=1)
		{stop("Parameter followup: c(0) No Fixed Follow-up, c(1,durfup) Fixed follow with duration durfup (durfup>0)")}
		
		### Fixed Follow-up 
		if (followup[1]==1 & length(followup)!=2)
		{stop("Parameter followup: c(0) No Fixed Follow-up, c(1,durfup) Fixed follow with duration durfup (durfup>0)")}
		
		### Fixed Follow-up and duration<=0
		if (followup[1]==1 & followup[2]<=0)
		{stop("Parameter followup: c(0) No Fixed Follow-up, c(1,durfup) Fixed follow with duration durfup (durfup>0)")}
		
		# Drop Out
		## drop out[1]== 0 | 1
		if (sum(dropout[1]==c(0,1))==0)
		{stop("Parameter Drop Out: c(0) No Drop out, c(1,gammae,gammac) drop out with hazard rate")}
		
		## not Drop Out and length>1
		if (dropout[1]==0 & length(dropout)!=1)
		{stop("Parameter Drop Out: c(0) No Drop out, c(1,gammae,gammac) drop out with hazard rate")}
		
		## Drop Out and length!=3
		if (dropout[1]==1 & length(dropout)!=3)
		{stop("Parameter Drop Out: c(0) No Drop out, c(1,gammae,gammac) drop out with hazard rate")}
		
		## Drop Out and Gamma<0
		if (dropout[1]==1 & sum(c(dropout[2]<0,dropout[3]<0))!=0)
		{stop("Parameter Drop Out: Hazard rate need to be positive")}
		
		# Duration
		## duration of study>0
		if (durstudy<=0)
		{stop("Duration of the Study (durstudy) need to be strictly positive")}
		
		## Accrual duration>0
		if (duraccrual<=0)
		{stop("Duration of the Accrual (duraccrual) need to be strictly positive")}
		
		# study duration>Acrual duration
		if (duraccrual>durstudy)
		{stop("Duration of the Study (durstudy) need to be greater than duration of the accrual (duraccrual)")}
		
		# Accrual accrual+Fixed Follow-up>study duration
		if (duraccrual+followup[2]-durstudy<0 & followup[1]==1)
		{stop("Duration of the accrual + Duration of the follow-up needs to pe greater or equal to duration of the study")}
		
		
	}
