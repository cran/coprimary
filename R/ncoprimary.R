############################################
#               ncoprimary.R               #
############################################


#'
#' @title Sample size calculation in clinical trials with two co-primary time-to-event endpoints 
#' @description 
#' ncoprimary() is used to calculate the sample size for phase III clinical trial with two co-primary endpoints to assess the efficacy of treatment between two groups.
#' @usage 
#' ncoprimary(design,Survhyp1,Survhyp2,alpha1,alpha2,duraccrual,durstudy,power,pe,
#' look,fup,dropout,dqol)
#' 
#' @param design Superiority=c(1,sided)[with sided=1 if 1-sided and 2 if 2-sided]; Non inferiority=c(2); Equivalence=c(3) 
#' @param Survhyp1 For Superiority=c(thyp,t,hype,Sc); for Non inferiority=c(thyp,t,hype,Sc,hypeA); for Equivalence=c(t,delta,Sc): parameters at time t for the first endpoint,
#' if thyp=1 then hype is survival rate or the rate of patients without HRQoL deterioration in experimental arm under the null hypothesis and hypeA is the survival rate in the experimental arm under the alternative hypothesis; 
#' if thyp=2 then hype is the hazard ratio under the null hypothesis and hypeA is the hazard ratio under the alternative hypothesis; Sc is survival rate in the control arm; 
#' delta is the log hazard ratio equivalence margin. When endpoint is HRQoL, the survival rate is replaced by the rate of patients without HRQoL deterioration.
#' @param Survhyp2 For Superiority=c(thyp,t,hype,Sc); for Non inferiority=c(thyp,t,hype,Sc,hypeA); for Equivalence=c(t,delta,Sc): parameters at time t for the second endpoint
#' if thyp=1 then hype is survival rate in experimental arm under the null hypothesis and hypeA is the survival rate in the experimental arm under the alternative hypothesis; 
#' if thyp=2 then hype is the hazard ratio under the null hypothesis and hypeA is the hazard ratio under the alternative hypothesis; Sc is survival rate in the control arm; 
#' delta is the log hazard ratio equivalence margin. When endpoint is HRQoL, the survival rate is replaced by the rate of patients without HRQoL deterioration.
#' @param alpha1 Type I error assigned to the first endpoint, for Non inferiority, Equivalence and 1-sided superiority is a vector of length one. 
#' For 2-sided superiority is a vector to length two c(alpha.low, alpha.up).
#' @param alpha2 Type I error assigned to the second endpoint, for Non inferiority, Equivalence and 1-sided superiority is a vector of length one. 
#' For 2-sided superiority is a vector to length two c(alpha.low, alpha.up).
#' @param duraccrual Accrual duration, expressed in number of days, months or years
#' @param durstudy Study duration, expressed in number of days, months or years
#' @param power 1-Probability of a type II error. Default value = 0.80.
#' @param pe Proportion (ratio) of patients assigned to the experimental arm (with 0<pe<1). Default value = 0.50.
#' @param look The number of interim analyses, c(1) for one final analysis; c(nb, bound, timing) for at least one interim analyses with bound=c(bound.eff,bound.fut):1-sided or bound=c(bound.low,bound.up):2-sided.
#' nb the number of planned looks, bound.eff and bound.fut corresponds to the type of boundaries used for efficacy (i.e. reject H0) and futility (i.e. reject H1). 
#' bound.fut=0: No futility monitoring, 1: Lan deMets O.Brien Fleming, 2: Lan deMets Pocock. 
#' bound.low and bound.up the type of lower and upper boundaries used (1: Lan deMets O.Brien Fleming, 2: Lan deMets Pocock, 3: O.Brien Fleming, 4: Pocock). Default value = 1.
#' @param fup Follow-up information, No fixed:c(0) (follow-up until the end of study); 
#' Fixed:c(1, durfollow) with durfollow is the duration of follow-up. Default value = 0.
#' @param dropout Drop out information, No drop out=c(0); Drop out=c(1,gammae,gammac) with gammae and gammac are the hazard drop out rates in experimental arm 
#' and control arm respectively. Default value = 0.
#' @param dqol number of targeted dimensions for the health related quality of life. Default value = 0.
#' 
#' @details The ncoprimary function computes the sample size for two primary endpoints. Both endpoints can be one time to event endpoint and health related quality of life (HRQoL) 
#' or two times to event endpoints.
#' 
#' @return Event: number of events estimated
#' @return Total: number of patients
#' @return Ne: number for experimental arm for each endpoint 
#' @return Nc: number for control arm for each endpoint 
#' @return HR: Hazard ratio for each endpoint
#' @examples 
#' 
#' ####################################################################################
#' ############ Design superiority:one-sided with two co-primary endpoints ############
#' ####################################################################################
#' ## - For endpoint 1: 3-year survival rates Se=0.75 and Sc=0.65, alpha1=0.02
#' ## - For endpoint 2: 4-year survival rates Se=0.70 and Sc=0.59, alpha2=0.03
#' ## with accrual duration of 2 years, study duration of 4 years and default values i.e 
#' ## power=0.80, pe=0.5, look=1, fup=0, dropout=0, dqol=0 
#' 
#' nc1 <- ncoprimary(design=c(1,1),Survhyp1=c(1,3,0.75,0.65),Survhyp2=c(1,4,0.70,0.59),
#' alpha1=0.02,alpha2=0.03,duraccrual=2,durstudy=4)
#' 
#' 
#' #####################################################################################
#'  ############ Design superiority:two-sided with two co-primary endpoints ############
#' #####################################################################################
#' ## - For endpoint 1: 2 target variables for the health related quality of life with 3-year 
#' ## rate without HRQoL deterioration Se=0.75 and Sc=0.67, alpha1=c(0.01,0.01)
#' ## - For endpoint 2: 4-year survival rates Se=0.86 and Sc=0.80, alpha2=c(0.015,0.015)
#' ## with accrual duration of 3 years, study duration of 6 years, power=0.90, look=c(2,c(1,1),0.5), 
#' ## and default values i.e  pe=0.5, fup=0, dropout=0
#' 
#' nc2 <- ncoprimary(design=c(1,2),Survhyp1=c(1,5,0.75,0.67),Survhyp2=c(1,5,0.86,0.80),
#' alpha1=c(0.01,0.01),alpha2=c(0.015,0.015),duraccrual=3,durstudy=6, power=0.90,
#' look=c(2,c(1,1),0.5),dqol=2)
#' 
#' 
#' #####################################################################################
#' ##############  Design non-inferiority with two co-primary endpoints ################
#' #####################################################################################
#' ## - For endpoint 1: 3-year survival rates Se=0.75 and Sc=SeA=0.75, alpha1=0.01
#' ## - For endpoint 2: 4-year survival rates Se=0.67 and Sc=SeA=0.80, alpha2=0.04
#' ## with accrual duration of 2 years, study duration of 6 years, power=0.95, pe=0.60 and 
#' ## default values i.e look=1, fup=0, dropout=0, dqol=0
#' 
#' nc3 <- ncoprimary(design=c(2),Survhyp1=c(1,4,0.65,0.75,0.75),Survhyp2=c(1,5,0.67,0.80,0.80),
#' alpha1=0.01,alpha2=0.04,duraccrual=2,durstudy=6,power=0.95,pe=0.60)
#' 
#' ####################################################################################
#' ################  Design superiority with two co-primary endpoints ################# 
#' ####################################################################################
#' 
#' ## - For endpoint 1: 2-year survival rate Sc=0.65 and log hazard equivalence margin delta=0.15 
#' ## and alpha1=0.025
#' ## - For endpoint 2: 1-year survival rate Sc=0.70 and log hazard equivalence margin delta=0.10 
#' ## and alpha2=0.025
#' ## with accrual duration of 3 years, study duration of 5 years, drop out hazard rate of 0.025 
#' ## per arm and default values i.e power=0.80, pe=0.5, look=1, fup=0, dqol=0 
#' 
#' nc4 <- ncoprimary(design=c(3),Survhyp1=c(2,0.15,0.65),Survhyp2=c(1,0.10,0.70),alpha1=0.025,
#' alpha2=0.025,duraccrual=3,durstudy=5,dropout=c(1,0.025,0.025))
#' 
#' 
#' @name ncoprimary
#' @export 

	
##'
##' @include datacheck.R
##' @include nbevent.R
##' @include probafix.R
##' @include probanofix.R
##' 
#NULL

 ncoprimary <- function(design,Survhyp1,Survhyp2,alpha1,alpha2,duraccrual,durstudy,power=0.80,pe=0.5,look=1,fup=0,dropout=0,dqol=0)

 {
	 alfa1=alpha1
	 alfa2=alpha2
	 beta=power
	 
	 datacheck(design=design,Survhyp=Survhyp1,pe=pe,alfa=alfa1,beta=beta,duraccrual=duraccrual,durstudy=durstudy,look=look,followup=fup,dropout=dropout)
	 #print("+ Parameters OK for endpoint 1")
	 datacheck(design=design,Survhyp=Survhyp2,pe=pe,alfa=alfa2,beta=beta,duraccrual=duraccrual,durstudy=durstudy,look=look,followup=fup,dropout=dropout)
	 #print("+ Parameters OK for endpoint 2")
	 
	 
	 if (design[1]==1 & design[2]==2)
	 {
		 alfa <- c(alfa1[1]+alfa2[1], alfa1[2]+alfa2[2])
	 }
	 else
	 {
		 alfa=alfa1+alfa2
	 }	
	 Sc1=nsurvival(design,Survhyp1,alpha1,duraccrual,durstudy,power,pe,look,fup,dropout,dqol)$Sc
	 Se1=nsurvival(design,Survhyp1,alpha1,duraccrual,durstudy,power,pe,look,fup,dropout,dqol)$Se
	 time1=nsurvival(design,Survhyp1,alpha1,duraccrual,durstudy,power,pe,look,fup,dropout,dqol)$time

	 
	 Sc2=nsurvival(design,Survhyp2,alpha2,duraccrual,durstudy,power,pe,look,fup,dropout,0)$Sc
	 Se2=nsurvival(design,Survhyp2,alpha2,duraccrual,durstudy,power,pe,look,fup,dropout,0)$Se
	 time2=nsurvival(design,Survhyp2,alpha2,duraccrual,durstudy,power,pe,look,fup,dropout,0)$time
	 
	 
	 print("########################################################")
	 print("|  SAMPLE SIZE CALCULATION FOR TWO SURVIVAL ENDPOINTS  |")
	 print("########################################################")
	 
	 
	 if (design[1]==1)
	 {
		 print("+-------------------------------------------------+")
		 print("| SURVIVAL SUPERIORITY TRIAL: TWO-SAMPLE TEST     |")
		 print("+-------------------------------------------------+")
	 }
	 if (design[1]==2)
	 {
		 print("+-------------------------------------------------+")
		 print("| SURVIVAL NON INFERIORITY TRIAL: TWO-SAMPLE TEST |")
		 print("+-------------------------------------------------+")
	 }
	 
	 if (design[1]==3)
	 {
		 print("+-------------------------------------------------+")
		 print("| SURVIVAL EQUIVALENCE TRIAL: TWO-SAMPLE TEST     |")
		 print("+-------------------------------------------------+")
	 }
	 print("+-----------------------------------------+")
	 print("|          TEST PARAMETERS                |")
	 print("+-----------------------------------------+")
	 if (design[1]==1 & design[2]==2)
	 {print("  - 1-Sided or 2-Sided Test: 2-Sided")}
	 else 
	 {print("  - 1-Sided or 2-Sided Test: 1-Sided")}
	 print(paste("  - Significance level alpha1 for endpoint 1:", round(sum(alfa1),3)))
	 print(paste("  - Significance level alpha2 for endpoint 2:",round(sum(alfa2),3)))
	 print(paste("  - Significance level alpha global        :",round(sum(alfa),3)))
	 
	 print(paste("  - Power:",power))
	 
	 print("+-----------------------------------------+")
	 print("|          STUDY PARAMETERS               |")
	 print("+-----------------------------------------+")
	 print(paste("  - Accrual Duration (duraccrual):",duraccrual))
	 
	 if (fup[1]==0)
	 {print("  - Follow-up: No Fixed Follow-up")}
	 else 
	 {print(paste("  - Follow-up: Fixed Follow-up of duration",fup[2]))}
	 
	 print(paste("  - Study Duration (durstudy):",durstudy))
	 print(paste("  - Assigned Fraction Experimental Arm:",pe))
	 
	 if (dropout[1]==0)
	 {print("  - Drop Out: No Drop Out")}
	 else
	 {
		 temp1 <- c("  - Drop Out: Drop Out with hazard rate (gammae,gammac)=(")
		 temp2 <- paste(dropout[2],",",dropout[3],")",sep="")
		 print(paste(temp1,temp2,sep=""))
	 }
	 
	 
	 print("+-----------------------------------------+")
	 print("|          INTERIM ANALYSIS               |")
	 print("+-----------------------------------------+")
	 print(paste("  - Number of Planned Analysis:",look[1]))
	 
	 if (look[1]!=1)
	 {
		 print(paste("  - Spacing of analysis:", paste(c(round(look[-c(1,2,3)],3),1),collapse=" ")))
		 boundary.txt <- c("Lan deMets O'Brien Fleming","Lan deMets Pocok","O'Brien Fleming","Pocok")
		if (design[1]==1)
		{
			if (design[2]==2)
			{
				print("  - Hypothesis to be rejected: Only H0")
				print(paste("  - Boundary to reject (H0- / H0+):",boundary.txt[look[2]],"/",boundary.txt[look[3]]),sep="")
			}
			if (design[2]==1 & look[3]==0)
			{
				print("  - Hypothesis to be rejected: Only H0")
				print(paste("  - Boundary to reject H0:",boundary.txt[look[2]]),sep="")
			}
			if (design[2]==1 & look[3]!=0)
			{
				print("  - Hypothesis to be rejected: H0 and H1")
				print(paste("  - Boundary to reject H0:",boundary.txt[look[2]]),sep="")
				print(paste("  - Boundary to reject H1:",boundary.txt[look[3]]),sep="")
			}
			if (design[2]==2 & alfa1[1]==alfa1[2] )
			{print(paste("  - Symmetric Boundary Alpha for endpoint 1: Lower=",alfa1[1],"/ Upper=",alfa1[2],sep=""))}
			if (design[2]==2 & alfa1[1]!=alfa1[2])
			{print(paste("  - Asymmetric Boundary Alpha for endpoint 1: Lower=",alfa1[1],"/ Upper=",alfa1[2],sep=""))}
			
			if (design[2]==2 & alfa2[1]==alfa2[2] )
			{print(paste("  - Symmetric Boundary Alpha for endpoint 2: Lower=",alfa2[1],"/ Upper=",alfa2[2],sep=""))}
			if (design[2]==2 & alfa2[1]!=alfa2[2])
			{print(paste("  - Asymmetric Boundary Alpha for endpoint 2: Lower=",alfa2[1],"/ Upper=",alfa2[2],sep=""))}
		}
	 }
	 
	 print("+-----------------------------------------+")
	 print("|    SURVIVAL PARAMETERS FOR ENDPOINT 1   |")
	 print("+-----------------------------------------+")
	 
	 print(paste("  - time (time) for endpoint 1:",time1))
	 print(paste("  - Survival Control (Sc) for endpoint 1:",Sc1))
	 
	 if(dqol!=0)
	 {
		 print(paste("  - Number of dimension quality of life:",dqol))
	 } 

		 if (design[1]==1)
		 {
			 HR1=nsurvival(design,Survhyp1,alpha1,duraccrual,durstudy,power,pe,look,fup,dropout,dqol)$HR
			 print(paste("  - Survival Experimental (Se) for endpoint 1:",round(Se1,3)))
			 print(paste("  - Hazard Ratio under Alternative Hypothesis (HR) for endpoint 1:",round(HR1,3)))
		 }
		 if (design[1]==2)
		 {
			 SeA1=nsurvival(design,Survhyp1,alpha1,duraccrual,durstudy,power,pe,look,fup,dropout,dqol)$SeA
			 HR1=nsurvival(design,Survhyp1,alpha1,duraccrual,durstudy,power,pe,look,fup,dropout,dqol)$HR
			 HR.alt1=nsurvival(design,Survhyp1,alpha1,duraccrual,durstudy,power,pe,look,fup,dropout,dqol)$alt.HR
			 
			 print(paste("  - Survival Non Inferiority Margin (Se1) for endpoint 1:",round(Se1,3)))
			 print(paste("  - Survival Experimental Arm Under Alternative Hypothesis (SeA) for endpoint 1:",round(SeA1,3)))
			 print(paste("  - Hazard Ratio under Null Hypothesis (HR) for endpoint 1:",round(HR1,2)))
			 print(paste("  - Hazard Ratio under Alternative Hypothesis (HR) for endpoint 1:",round(HR.alt1,3)))
		 }
		 
		 if (design[1]==3)
		 {
			 logHR1=nsurvival(design,Survhyp1,alpha1,duraccrual,durstudy,power,pe,look,fup,dropout,dqol)$logHR
			 
			 print(paste("  - Equivalence Margin log(Hazard Ratio) for endpoint 1:",logHR1))
			 print(paste("  - One-sided hypothesis 1 for endpoint 1: H01: h>=",round(exp(Survhyp1[3]),3)," vs H11: h<",round(exp(Survhyp1[3]),3),sep=""))
			 print(paste("  - One-sided hypothesis 2 for endpoint 1: H02: h<=",round(exp(-Survhyp1[3]),3)," vs H12: h>",round(exp(-Survhyp1[3]),3),sep=""))
		 }
		 events1=nsurvival(design,Survhyp1,alpha1,duraccrual,durstudy,power,pe,look,fup,dropout,dqol)$events
		 n1=nsurvival(design,Survhyp1,alpha1,duraccrual,durstudy,power,pe,look,fup,dropout,dqol)$subjects
		 interim1=nsurvival(design,Survhyp1,alpha1,duraccrual,durstudy,power,pe,look,fup,dropout,dqol)$interim

	 
	 
	 print("+-----------------------------------------+")
	 print("|    SURVIVAL PARAMETERS FOR ENDPOINT 2   |")
	 print("+-----------------------------------------+")
	 print(paste("  - time (time) for endpoint 2:",time2))
	 print(paste("  - Survival Control (Sc) for endpoint 2:",Sc2))
	 if (design[1]==1)
	 {
		 HR2=nsurvival(design,Survhyp2,alpha2,duraccrual,durstudy,power,pe,look,fup,dropout,0)$HR
		 print(paste("  - Survival Experimental (Se) for endpoint 2:",round(Se2,3)))
		 print(paste("  - Hazard Ratio under Alternative Hypothesis (HR) for endpoint 2:",round(HR2,3)))
	 }
	 if (design[1]==2)
	 {
		 SeA2=nsurvival(design,Survhyp2,alpha2,duraccrual,durstudy,power,pe,look,fup,dropout,0)$SeA
		 HR2=nsurvival(design,Survhyp2,alpha2,duraccrual,durstudy,power,pe,look,fup,dropout,0)$HR
		 HR.alt2=nsurvival(design,Survhyp2,alpha2,duraccrual,durstudy,power,pe,look,fup,dropout,0)$alt.HR
		 
		 print(paste("  - Survival Non Inferiority Margin (Se) for endpoint 2:",round(Se2,3)))
		 print(paste("  - Survival Experimental Arm Under Alternative Hypothesis (SeA) for endpoint 2:",round(SeA2,3)))
		 print(paste("  - Hazard Ratio under Null Hypothesis (HR) for endpoint 2:",round(HR2,2)))
		 print(paste("  - Hazard Ratio under Alternative Hypothesis (HR) for endpoint 2:",round(HR.alt2,3)))
	 }
	 
	 if (design[1]==3)
	 {
		 logHR2=nsurvival(design,Survhyp2,alpha2,duraccrual,durstudy,power,pe,look,fup,dropout,0)$logHR
		 
		 print(paste("  - Equivalence Margin log(Hazard Ratio) for endpoint 2:",logHR2))
		 print(paste("  - One-sided hypothesis 1 for endpoint 2: H01: h>=",round(exp(Survhyp2[3]),3)," vs H11: h<",round(exp(Survhyp2[3]),3),sep=""))
		 print(paste("  - One-sided hypothesis 2 for endpoint 2: H02: h<=",round(exp(-Survhyp2[3]),3)," vs H12: h>",round(exp(-Survhyp2[3]),3),sep=""))
	 }
	 
	  
	 events2=nsurvival(design,Survhyp2,alpha2,duraccrual,durstudy,power,pe,look,fup,dropout,0)$events
	 n2=nsurvival(design,Survhyp2,alpha2,duraccrual,durstudy,power,pe,look,fup,dropout,0)$subjects
	 interim2=nsurvival(design,Survhyp2,alpha2,duraccrual,durstudy,power,pe,look,fup,dropout,0)$interim
	 
	 if(n1>=n2)
	 {
		 events=events1
		 n=n1
	 	interim=interim1
	 }
	 else{ 
		 events=events2
		 n=n2 
		 interim=interim2
	 }
	 
	 
	 print("+-----------------------------------------+")
	 print("|          SAMPLE SIZE                    |")
	 print("+-----------------------------------------+")
	 
	 print(paste("  - Number of Events:",round(events)))
	 print(paste("  - Total Number of Subjects:",round(n)))
	 print(paste("  - Number of Subjects (Arm Exp., Arm Contr.): (",round(n*pe),",",round(n*(1-pe)),")",sep=""))
	 
	if (look[1]>1)
	{
		print("+ Boundaries ")
		print(interim)

	}
	 
	 list(alpha=alfa,power=power,accrual=duraccrual,look=look,followup=fup,durstuty=durstudy, events=events,totalnumber=n, Ne=n/2, Nc=n/2)
	 
 }	