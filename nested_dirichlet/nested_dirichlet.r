library(MASS)
library(mvtnorm) 
library(LaplacesDemon)
library(matrixStats)

set.seed(12345)

# rm(list = ls())
# #
# #########################
# #########################
# # simulate the data set from different normal distributions
# #########################
# #########################
# #
# rDiscreta<-function(p){
#  u<-runif(1)
#  P<-cumsum(p)
#  val<-sum(P<u)+1
#  return(val)}
# #
# K<-4 # number of different distributions (clusters)
# J<-40 # number of elements in the sample
# Ij<-rep(100,J) # number of replication for each element
# #
# meanvec1<-c(-3,3)
# sigmavec1<-c(2,2)
# p1<-c(0.30,0.70)
# meanvec2<-c(3)
# sigmavec2<-c(9)
# p2<-1
# meanvec3<-c(-3)
# sigmavec3<-c(9)
# p3<-1
# meanvec4<-c(-5,5,10)
# sigmavec4<-c(2,2,2)
# p4<-c(0.25,0.65,0.10)
# #
# id<-numeric(sum(Ij)) # specifies from each distribution each element comes from
# iddentr<-numeric(sum(Ij))
# Y<-numeric(sum(Ij)) # observations
# #
# set.seed(1034)
# cont<-1
# for (j in 1:40){
# 	dist<-rDiscreta(rep(1/K,K))
# 	for (i in 1:Ij[j]){
# 		id[cont]<-dist
# 		if (id[cont]==1){
# 			iddentr[cont]<-rDiscreta(p1)
# 			Y[cont]<-rnorm(1,meanvec1[iddentr[cont]],sqrt(sigmavec1[iddentr[cont]]))}
# 		if (id[cont]==2){
# 			iddentr[cont]<-rDiscreta(p2)
# 			Y[cont]<-rnorm(1,meanvec2[iddentr[cont]],sqrt(sigmavec2[iddentr[cont]]))}	
# 		if (id[cont]==3){
# 			iddentr[cont]<-rDiscreta(p3)
# 			Y[cont]<-rnorm(1,meanvec3[iddentr[cont]],sqrt(sigmavec3[iddentr[cont]]))}		
# 		if (id[cont]==4){
# 			iddentr[cont]<-rDiscreta(p4)
# 			Y[cont]<-rnorm(1,meanvec4[iddentr[cont]],sqrt(sigmavec4[iddentr[cont]]))}
# 		cont<-cont+1}}
# #
# tauj_verd<-id # specifies from each distribution each element comes from in level of elements and their replication
# Sj_verd<-numeric(J) # specifies from each distribution each element comes from in level of elements
# for (i in 1:J){
# 	Sj_verd[i]<-as.numeric(as.character(data.frame(table(tauj_verd[(cumsum(Ij)[i]-Ij[i]+1):cumsum(Ij)[i]]))[,1]))}
# #
# id<-rep(1,Ij[1])
# for (j in 2:J) id<-c(id,rep(j,Ij[j]))		
# dataset<-matrix(c(id,Y),ncol=2) # final data set
#
#########################
#########################
# useful functions
#########################
#########################
#
############
# sample from a Normal Inverse gamma distribituion (mu, lambda, alpha, beta)
############
rinvgamma<-function(media.mi,lambda,alpha,beta){
	n<-length(media.mi)
	sigma2<-1/(rgamma(n,alpha,beta))
	mu<-rnorm(n,media.mi,sqrt(sigma2/lambda))
	sample<-cbind(mu,sigma2)
	return(sample)}
#
############
# sample from a Multinomial distribituion (p)
############
rDiscreta<-function(p){
 u<-runif(1)
 P<-cumsum(p)
 val<-sum(P<u)+1
 return(val)}
#
############
# calculate the density function of a student t distribituion (yij, ni, mu, sigma)
############
dstudentt<-function(yij,ni,mu,desvio){
 dens<-(gamma((ni+1)/2)/gamma(ni/2))*((1+(1/ni)*((yij-mu)/desvio)**2)**(-(ni+1)/2))*(1/(sqrt(pi*ni)*desvio))
 return(dens)}
#
############
# calculate the marginalized likelihood of Yj given (Xij ,tauj) and specific l  -> x_ij=mu_ij,  tau_j=sigma2_j
############
dmarglikeli<-function(Yjl,lambda,agam,bgam,mu){
 mlkj<-length(Yjl)
 media<-mean(Yjl)
 varia<-var(Yjl)*(mlkj-1)
 if (is.na(varia)) varia<-0
 logdens<-(lgamma(mlkj/2+agam)-lgamma(agam))+((-mlkj/2)*log(2*pi))+(0.5*(log(lambda)-log(mlkj+lambda)))+(agam*log(bgam))+((-(mlkj/2+agam))*log(bgam+(varia/2)+((lambda*mlkj*((media-mu)**2))/(2*(mlkj+lambda)))))
 return(logdens)} 
#
############
# sample from conditional a posteriori distribution of theta for all k and l
############	
posteriori.theta<-function(dataset, Zij, tauj, lambda, agam, bgam,mu, Xij, b){
	mlk<-aggregate(rep(1,nrow(dataset)), by = list(Zij, tauj), FUN = "sum")[,3]
#	Ymeanlk<-aggregate(dataset[,2], by = list(Zij, tauj), FUN = "mean")[,3]
#	Yvarlk<-aggregate(dataset[,2], by = list(Zij, tauj), FUN = "var")[,3]*(mlk-1)
	dataset2=dataset
	dataset2$scarto=rep(0,dim(dataset)[1])
	for (r in 1:dim(dataset)[1]){
	dataset2$scarto[r] = t(Xij[r,])%*%b[dataset[r,1],]
	}
	dataset2$PV1MATH=dataset2$PV1MATH-dataset2$scarto
	Ymeanlk<-aggregate(dataset2[,2], by = list(Zij, tauj), FUN = "mean")[,3]
	Yvarlk<-aggregate(dataset2[,2], by = list(Zij, tauj), FUN = "var")[,3]*(mlk-1)
	mupost<-((mlk*Ymeanlk)+(lambda*mu))/(lambda+mlk)
	lambdapost<-lambda+mlk
	agampost<-agam+(mlk/2)
	Yvarlk[is.na(Yvarlk)] <- 0
	bgampost<-bgam+(Yvarlk+((mlk*lambda*(Ymeanlk-mu)**2)/(lambda+mlk)))/2
	theta.samp<-rinvgamma(mupost,lambdapost,agampost,bgampost)
	return(theta.samp)}
#--------------------
############
# sample from conditional a posteriori distribution of Zij =rij
############	
posteriori.Zij<-function(dataset, Zij, tauj, K, theta.samp, ni, mu, desvio, Lk, beta,Xij, b){
	for (k in 1:K){
		observ<-which(tauj==k)
		dataset2=dataset
		dataset2$scarto=rep(0,dim(dataset)[1])
		for (r in 1:dim(dataset)[1]){
		  dataset2$scarto[r] = t(Xij[r,])%*%b[dataset[r,1],]
		}
		Yclusterk<-dataset[observ,2]-dataset2$scarto[observ]
		Zclusterk<-Zij[observ]
		cont<-1
		for (i in observ){
			probZij<-numeric(Lk[k])
			contmlk<-Zclusterk[-cont]
			l<-seq(1:Lk[k])
			mlk<-rep(0,Lk[k])
			need<-data.frame(table(contmlk[which(contmlk<=Lk[k])]))
			#modifico:
			#if (ncol(need) >= 2) { 
			
			categr <- as.numeric(as.character(need[,1]))
			if (nrow(need) != 0) mlk[categr] <- as.numeric(need[,2])
			#} else {
			#  mlk <- rep(0, Lk[Sj[j]])  # Inizializza mlk per evitare errore
			#}
			#
			#categr<-as.numeric(as.character(need[,1]))
			#mlk[categr]<-as.numeric(need[,2])
			probZij<-mlk*dnorm(Yclusterk[cont],theta.samp[(cumsum(Lk)[k]-Lk[k]+l),1],sqrt(theta.samp[(cumsum(Lk)[k]-Lk[k]+l),2]))
			probZij<-c(probZij,beta*dstudentt(Yclusterk[cont],ni,mu,desvio))		
			probZij<-probZij/sum(probZij)
			Zij[i]<-rDiscreta(probZij)
			Zclusterk<-Zij[observ]
			#print(Zclusterk)
			cont<-cont+1}
		while (length(table(Zclusterk))<max(Zclusterk)){ # exclude empty clusters
			categr<-as.numeric(as.character(data.frame(table(Zclusterk))[,1]))
			categd<-seq(1:length(table(Zclusterk)))
			dif<-which(categr!=categd)
			for (i in 1:length(Zclusterk)) if (Zclusterk[i]>dif[1]) Zclusterk[i]<-Zclusterk[i]-1}
		Zij[observ]<-Zclusterk}
	return(Zij)}
#
############
# sample from conditional a posteriori distribution of Zij for just one j
############	
posteriori.Zij2<-function(dataset, Zij, tauj, K, theta.samp, ni, mu, desvio, Lk, beta, j, Sj,Xij,b){
	k<-Sj[j]
	observ<-which(tauj==k)
	dataset2=dataset
	dataset2$scarto=rep(0,dim(dataset)[1])
	for (r in 1:dim(dataset)[1]){
	  dataset2$scarto[r] = t(Xij[r,])%*%b[dataset[r,1],]
	}
	Yclusterk<-dataset[observ,2]-dataset2$scarto[observ]
	Zclusterk<-Zij[observ]
	cont<-1
	for (i in observ){
		probZij<-numeric(Lk[k])
		contmlk<-Zclusterk[-cont]
		l<-seq(1:Lk[k])
		mlk<-rep(0,Lk[k])
		need<-data.frame(table(contmlk[which(contmlk<=Lk[k])]))
		#modifico:
		#if (ncol(need) >= 2) { 
		categr <- as.numeric(as.character(need[,1]))
		if (nrow(need) != 0) mlk[categr] <- as.numeric(need[,2])
		#} else {
		 # mlk <- rep(0, Lk[Sj[j]])  # Inizializza mlk per evitare errore
		#}
		#
		#categr<-as.numeric(as.character(need[,1]))
		#mlk[categr]<-as.numeric(need[,2])
		probZij<-log(mlk)+dnorm(Yclusterk[cont],theta.samp[(cumsum(Lk)[k]-Lk[k]+l),1],sqrt(theta.samp[(cumsum(Lk)[k]-Lk[k]+l),2]),log=TRUE)
		probZij<-c(probZij,log(beta)+dst(Yclusterk[cont],mu=mu,sigma=desvio,nu=ni,log=TRUE))	
		#print(Yclusterk[cont])
		probZij<-exp(probZij-logSumExp(probZij))
		Zij[i]<-rDiscreta(probZij)
		Zclusterk<-Zij[observ]
		cont<-cont+1}
	while (length(table(Zclusterk))<max(Zclusterk)){ # exclude empty clusters
		categr<-as.numeric(as.character(data.frame(table(Zclusterk))[,1]))
		categd<-seq(1:length(table(Zclusterk)))
		dif<-which(categr!=categd)
		for (i in 1:length(Zclusterk)) if (Zclusterk[i]>dif[1]) Zclusterk[i]<-Zclusterk[i]-1}
	Zij[observ]<-Zclusterk
	return(Zij)}
#
############
# build a proposal of new cluster for each Sj and realocate its Zij---------------
############
sample.tau.xi<-function(dataset, Ij, Sj, Zij, tauj, K, Lk, theta.samp, j, alpha, beta, ni, mu, desvio,Xij,b){
	nk<-numeric(K)
	for (k in 1:K) nk[k]<-sum(Sj[-j]==k) #numerosità ogni cluster senza paese j
	observ<-which(dataset[,1]==j)
	dataset2=dataset
	dataset2$scarto=rep(0,dim(dataset)[1])
	for (r in 1:dim(dataset)[1]){
	  dataset2$scarto[r] = t(Xij[r,])%*%b[dataset[r,1],]
	}
	Yclusterj<-dataset[observ,2]-dataset2$scarto
	Zclusterj<-Zij[observ]
	nk[Sj[j]]<-0 # we force Sj candidate to be different of the current value
	nk<-c(nk,alpha)
	probSj<-nk/sum(nk)
	Sprop<-rDiscreta(probSj)
	Sold<-Sj[j]
	Zprop<-numeric(Ij[j])
	Zold<-numeric(Ij[j])
	lprior<-0
	lfunctrans<-0
	lpriorold<-0
	lfunctransold<-0
	taujcand<-tauj
	taujcand[(cumsum(Ij)[j]-Ij[j]+1):cumsum(Ij)[j]]<-rep(K+1,Ij[j])
	nj<-sum(Sj[-j]==Sold)
	if (Sprop <= K){
		for (i in 1:Ij[j]){
			l<-1:Lk[Sprop]
			Zvector<-c(Zij[which(tauj==Sprop)],Zprop[which(Zprop<=Lk[Sprop] & Zprop>0)])
			mlk<-table(Zvector)
			probZij<-mlk*dnorm(Yclusterj[i],theta.samp[(cumsum(Lk)[Sprop]-Lk[Sprop]+l),1],sqrt(theta.samp[(cumsum(Lk)[Sprop]-Lk[Sprop]+l),2]))
			if (length(which(Zprop>Lk[Sprop]))>0){
				mlk2<-table(Zprop[which(Zprop>Lk[Sprop])])
				probZij2<-mlk2*dstudentt(Yclusterj[i],ni,mu,desvio)
				mlk<-c(mlk,mlk2)
				probZij<-c(probZij,probZij2)}
			mlk<-c(mlk,beta)
			probZij<-c(probZij,mlk[length(mlk)]*dstudentt(Yclusterj[i],ni,mu,desvio))		
			probZij<-probZij/sum(probZij)
			mlk<-mlk/sum(mlk)
			Zprop[i]<-rDiscreta(probZij)
			lprior<-lprior+log(mlk[Zprop[i]])
			lfunctrans<-lfunctrans+log(probZij[Zprop[i]])
			if (nj > 0){
				l<-1:Lk[Sold]
				Zvector<-c(Zij[which(taujcand==Sold)],Zold[which(Zold<=Lk[Sold] & Zold>0)])
				mlk<-rep(0,Lk[Sold])
				need<-data.frame(table(Zvector))
				#modifico:
				#if (ncol(need) >= 2) { 
				categr <- as.numeric(as.character(need[,1]))
				if (nrow(need) != 0) mlk[categr] <- as.numeric(need[,2])
				#} else {
				#  mlk <- rep(0, Lk[Sj[j]])  # Inizializza mlk per evitare errore
				#}
				#
				#categr<-as.numeric(as.character(need[,1]))
				#mlk[categr]<-as.numeric(need[,2])
				probZij<-mlk*dnorm(Yclusterj[i],theta.samp[(cumsum(Lk)[Sold]-Lk[Sold]+l),1],sqrt(theta.samp[(cumsum(Lk)[Sold]-Lk[Sold]+l),2]))
				if (length(which(mlk==0))==0){
					mlk<-c(mlk,beta)
					probZij<-c(probZij,mlk[length(mlk)]*dstudentt(Yclusterj[i],ni,mu,desvio))} else {
					empty<-which(mlk==0)
					#print(paste("1:", mlk, Zclusterj[i]))
					if (length(mlk) >= Zclusterj[i] && mlk[Zclusterj[i]]==0){
						mlk[Zclusterj[i]]<-beta
						probZij[Zclusterj[i]]<-mlk[Zclusterj[i]]*dstudentt(Yclusterj[i],ni,mu,desvio)} else {
						mlk[empty[1]]<-beta
						probZij[empty[1]]<-mlk[empty[1]]*dstudentt(Yclusterj[i],ni,mu,desvio)}}	
				probZij<-probZij/sum(probZij)
				mlk<-mlk/sum(mlk)
				Zold[i]<-Zclusterj[i]
				lpriorold<-lpriorold+log(mlk[Zclusterj[i]])
				lfunctransold<-lfunctransold+log(probZij[Zclusterj[i]])} else {
				Zold[1]<-Zclusterj[1]
				Zvector<-c(Zij[which(taujcand==Sold)],Zold[which(Zold <= Lk[Sold] & Zold>0)])
				mlk<-rep(0,Lk[Sold])
				need<-data.frame(table(Zvector))
			  categr<-as.numeric(as.character(need[,1]))
				if(length(categr)>0 && ncol(need)>2) mlk[categr]<-as.numeric(need[,2])
				probZij<-mlk*dstudentt(Yclusterj[i],ni,mu,desvio)
				if (length(which(mlk==0))==0){
					mlk<-c(mlk,beta)
					probZij<-c(probZij,mlk[length(mlk)]*dstudentt(Yclusterj[i],ni,mu,desvio))}	else { 
					empty<-which(mlk==0)
					#print(paste("2:", mlk, Zclusterj[i]))
					
					if (length(mlk) >= Zclusterj[i] && mlk[Zclusterj[i]]==0){ 
						mlk[Zclusterj[i]]<-beta
						probZij[Zclusterj[i]]<-mlk[Zclusterj[i]]*dstudentt(Yclusterj[i],ni,mu,desvio)} else {
						mlk[empty[1]]<-beta
						probZij[empty[1]]<-mlk[empty[1]]*dstudentt(Yclusterj[i],ni,mu,desvio)}}	
				probZij<-probZij/sum(probZij)				
				mlk<-mlk/sum(mlk)
				Zold[i]<-Zclusterj[i]
				lpriorold<-lpriorold+log(mlk[Zclusterj[i]])
				lfunctransold<-lfunctransold+log(probZij[Zclusterj[i]])}}}				
	if (Sprop > K){
		Zprop[1]<-1
		for (i in 2:Ij[j]){
			mlk<-table(Zprop[which(Zprop>0)])
			probZij<-mlk*dstudentt(Yclusterj[i],ni,mu,desvio)
			mlk<-c(mlk,beta)
			probZij<-c(probZij,mlk[length(mlk)]*dstudentt(Yclusterj[i],ni,mu,desvio))		
			probZij<-probZij/sum(probZij)
			mlk<-mlk/sum(mlk)
			Zprop[i]<-rDiscreta(probZij)
			lprior<-lprior+log(mlk[Zprop[i]])
			lfunctrans<-lfunctrans+log(probZij[Zprop[i]])
			if (nj > 0){
				l<-1:Lk[Sold]
				Zvector<-c(Zij[which(taujcand==Sold)],Zold[which(Zold<=Lk[Sold] & Zold>0)])
				mlk<-rep(0,Lk[Sold])
				need<-data.frame(table(Zvector))
				#mofidico:
				#if (ncol(need) >= 2) {
				categr <- as.numeric(as.character(need[,1]))
				if (nrow(need) != 0) mlk[categr] <- as.numeric(need[,2])
				#} else {
				 # mlk <- rep(0, Lk[Sj[j]])  # Inizializza mlk per evitare errore
				#}
				#
				#categr<-as.numeric(as.character(need[,1]))
				#mlk[categr]<-as.numeric(need[,2])
				probZij<-mlk*dnorm(Yclusterj[i],theta.samp[(cumsum(Lk)[Sold]-Lk[Sold]+l),1],sqrt(theta.samp[(cumsum(Lk)[Sold]-Lk[Sold]+l),2]))
				if (length(which(mlk==0))==0){
					mlk<-c(mlk,beta)
					probZij<-c(probZij,mlk[length(mlk)]*dstudentt(Yclusterj[i],ni,mu,desvio))} else {
					empty<-which(mlk==0)
					#print(paste("3:", mlk, Zclusterj[i]))
					
					if (length(mlk) >= Zclusterj[i] && mlk[Zclusterj[i]]==0){
						mlk[Zclusterj[i]]<-beta
						probZij[Zclusterj[i]]<-mlk[Zclusterj[i]]*dstudentt(Yclusterj[i],ni,mu,desvio)} else {
						mlk[empty[1]]<-beta
						probZij[empty[1]]<-mlk[empty[1]]*dstudentt(Yclusterj[i],ni,mu,desvio)}}	
				probZij<-probZij/sum(probZij)
				mlk<-mlk/sum(mlk)
				Zold[i]<-Zclusterj[i]
				lpriorold<-lpriorold+log(mlk[Zclusterj[i]])
				lfunctransold<-lfunctransold+log(probZij[Zclusterj[i]])} else {
				Zold[1]<-Zclusterj[1]
				Zvector<-c(Zij[which(taujcand==Sold)],Zold[which(Zold <= Lk[Sold] & Zold>0)])
				mlk<-rep(0,Lk[Sold])
				need<-data.frame(table(Zvector))
				#modifichiamo le prox due righe:
				#if (ncol(need) >= 2) {  
				categr <- as.numeric(as.character(need[,1]))
				if (nrow(need) != 0) mlk[categr] <- as.numeric(need[,2])
				#} else {
				#  mlk <- rep(0, Lk[Sold])  # Inizializza mlk per evitare errore
				#}
				#
				#categr<-as.numeric(as.character(need[,1]))
				if(length(categr)>0 && ncol(need)>2) mlk[categr]<-as.numeric(need[,2])
				probZij<-mlk*dstudentt(Yclusterj[i],ni,mu,desvio)
				if (length(which(mlk==0))==0){
					mlk<-c(mlk,beta)
					probZij<-c(probZij,mlk[length(mlk)]*dstudentt(Yclusterj[i],ni,mu,desvio))}	else { 
					empty<-which(mlk==0)
					#print(paste("4:", mlk, Zclusterj[i]))
					
					if (length(mlk) >= Zclusterj[i] && mlk[Zclusterj[i]]==0){
						mlk[Zclusterj[i]]<-beta
						probZij[Zclusterj[i]]<-mlk[Zclusterj[i]]*dstudentt(Yclusterj[i],ni,mu,desvio)} else {
						mlk[empty[1]]<-beta
						probZij[empty[1]]<-mlk[empty[1]]*dstudentt(Yclusterj[i],ni,mu,desvio)}}	
				probZij<-probZij/sum(probZij)				
				mlk<-mlk/sum(mlk)
				Zold[i]<-Zclusterj[i]
				lpriorold<-lpriorold+log(mlk[Zclusterj[i]])
				lfunctransold<-lfunctransold+log(probZij[Zclusterj[i]])}}}
	Sjcand<-Sj
	Sjcand[j]<-Sprop
	taujcand<-tauj
	taujcand[(cumsum(Ij)[j]-Ij[j]+1):cumsum(Ij)[j]]<-rep(Sprop,Ij[j])
	Zijcand<-Zij
	Zijcand[(cumsum(Ij)[j]-Ij[j]+1):cumsum(Ij)[j]]<-Zprop
	return(list(Sjcand, taujcand, Zijcand, lprior, lfunctrans, Sprop, Zprop, Yclusterj, Zclusterj,lpriorold,lfunctransold))}
#
############
# calculate the acceptance rate of the pair Sj and all Zij
############	
prob.accept<-function(dataset, Yclusterj, Sj, Zij, tauj, Lk, theta.samp, j, Zijcand, taujcand, Sprop, Zprop, lambda, mu, agam, bgam, lprior, lfunctrans, lpriorold, lfunctransold,Xij,b){
	marglikel<-0
	marglikelold<-0
	#
	dataset2=dataset
	dataset2$scarto=rep(0,dim(dataset)[1])
	for (r in 1:dim(dataset)[1]){
	  dataset2$scarto[r] = t(Xij[r,])%*%b[dataset[r,1],]
	}
	dataset2$PV1MATH=dataset2$PV1MATH-dataset2$scarto
	Yclusterk<-dataset[which(tauj==Sprop),2]-dataset2$scarto[which(tauj==Sprop)]
	Zclusterk<-Zij[which(tauj==Sprop)]
	Lsprop<-length(table(c(Zclusterk,Zprop)))
	if (length(Zclusterk)==0){
		mlk<-rep(0,Lsprop)
		Ymeanlk<-rep(0,Lsprop)
		Yvarlk<-rep(0,Lsprop)}
	if (length(Zclusterk)>0){
		mlk<-table(Zclusterk)
		Ymeanlk<-aggregate(Yclusterk, by = list(Zclusterk,rep(1,length(Zclusterk))), FUN = "mean")[,3]
		Ymeanlk[is.na(Ymeanlk)] <- 0
		Yvarlk<-aggregate(Yclusterk, by = list(Zclusterk,rep(1,length(Zclusterk))), FUN = "var")[,3]*(mlk-1)
		Yvarlk[is.na(Yvarlk)] <- 0
		dif<-Lsprop-length(mlk)
		if (dif > 0){
			mlk<-c(mlk,rep(0,dif))
			Ymeanlk<-c(Ymeanlk,rep(0,dif))
			Yvarlk<-c(Yvarlk,rep(0,dif))}}
	mupost<-((mlk*Ymeanlk)+(lambda*mu))/(lambda+mlk)
	lambdapost<-lambda+mlk
	agampost<-agam+(mlk/2)
	bgampost<-bgam+(Yvarlk+((mlk*lambda*(Ymeanlk-mu)**2)/(lambda+mlk)))/2
	for (l in 1:Lsprop){
		Yclusterjl<-Yclusterj[Zprop==l]
		if (length(Yclusterjl)>0) marglikel<-marglikel+dmarglikeli(Yclusterjl,lambdapost[l],agampost[l],bgampost[l],mupost[l])}
	#
	Yclusterk<-dataset2[which(taujcand==Sj[j]),2]
	Zclusterk<-Zijcand[which(taujcand==Sj[j])]
	if (length(Zclusterk)==0){
		mlk<-rep(0,Lk[Sj[j]])
		Ymeanlk<-rep(0,Lk[Sj[j]])
		Yvarlk<-rep(0,Lk[Sj[j]])}
	if (length(Zclusterk)>0){
		mlk<-table(Zclusterk)
		Ymeanlk<-aggregate(Yclusterk, by = list(Zclusterk,rep(1,length(Zclusterk))), FUN = "mean")[,3]
		Ymeanlk[is.na(Ymeanlk)] <- 0
		Yvarlk<-aggregate(Yclusterk, by = list(Zclusterk,rep(1,length(Zclusterk))), FUN = "var")[,3]*(mlk-1)
		Yvarlk[is.na(Yvarlk)] <- 0
		dif<-Lk[Sj[j]]-length(mlk)
		if (dif > 0){
			mlk<-c(mlk,rep(0,dif))
			Ymeanlk<-c(Ymeanlk,rep(0,dif))
			Yvarlk<-c(Yvarlk,rep(0,dif))}}
	mupost<-((mlk*Ymeanlk)+(lambda*mu))/(lambda+mlk)
	lambdapost<-lambda+mlk
	agampost<-agam+(mlk/2)
	bgampost<-bgam+(Yvarlk+((mlk*lambda*(Ymeanlk-mu)**2)/(lambda+mlk)))/2
	Zold<-Zij[(cumsum(Ij)[j]-Ij[j]+1):cumsum(Ij)[j]]
	for (l in 1:Lk[Sj[j]]){
		Yclusterjl<-Yclusterj[Zold==l]
		if (length(Yclusterjl)>0) marglikelold<-marglikelold+dmarglikeli(Yclusterjl,lambdapost[l],agampost[l],bgampost[l],mupost[l])}	
	#
	paccept<-exp(marglikel+lprior+lfunctransold-(marglikelold+lpriorold+lfunctrans))
	#modifichiamo noi
	if (is.na(paccept)){
	  paccept=1
	  print(marglikel+lprior+lfunctransold-(marglikelold+lpriorold+lfunctrans))
	}
	return(list(paccept,lprior,lpriorold,lfunctrans,lfunctransold,marglikel,marglikelold))}
#
#
#########################
#########################
# Model initialization
#########################
#########################
#
########### data set needs to be ordered by elements and replications of each element
########### colunm 1 identifies the elements and the colunm 2 has the observations for each element
#
#setwd("/Users/gaiacaringi/Documents/Project\ Bayesian")
caminho<-"/Users/gaiacaringi/Documents/Project\ Bayesian" # this directory is where the files with MCMC output will be saved, then you need to change it for a directory of your computer
#data=read.csv("OECD_data.csv")
#data = data %>% group_by(CNT) %>% filter(n()>3) %>% ungroup()
#data = as.data.frame(data)
#data$CNT <- as.numeric(as.factor(data$CNT))
dataset <- data[,c(2,13)]
Xij=data[,-c(1,2,3,7,9,10,11,12,13)]
Xij=as.matrix(Xij)
J<-max(dataset[,1])
Nt<-nrow(dataset)
alpha<-1 # total mass parameter
beta<-1 # total mass parameter
agam<-3 # alpha hyperparameter of Inverse gamma distribution (G0)
bgam<-5 # beta hyperparameter of Inverse gamma distribution (G0)
mu<-0 # mu hyperparameter of Inverse gamma distribution (G0)
lambda<-0.01 # lambda hyperparameter of Inverse gamma distribution (G0)
ni<-2*agam # parameter of student-t distribution
desvio<-sqrt((2*bgam*(1+lambda))/(ni*lambda)) # parameter of student-t distribution
Ij<-numeric(J)
for (i in 1:J) Ij[i]<-sum(dataset[,1]==i)
#
# Initialize Sj, Zij and theta, beta
#
library(dplyr)
sigmab <- 1
bvec_totale <- matrix(0, nrow = J, ncol = dim(Xij)[2])

for (j in 1:J) {
  d <- sum(dataset$CNT == j)  # Conta quanti elementi hanno CNT == j
  bvec <- rnorm(dim(Xij)[2], 0, sigmab)  # Genera d numeri casuali
  bvec_totale[j,] <- bvec  # Assegna i valori ai posti giusti
}
b=bvec_totale

Sj<-rep(1,J) # distributional cluster membership indicator of elements
Zij<-rep(1,Nt) # observational cluster membership indicator for the replications of elements
K<-length(table(Sj))
Lk<-numeric(K)
tauj<-numeric()
for (j in 1:length(Sj)) tauj<-c(tauj,rep(Sj[j],Ij[j]))
for (k in 1:K) Lk[k]<-length(table(Zij[which(tauj==k)]))
set.seed(1000)
theta.samp<-posteriori.theta(dataset, Zij, tauj, lambda, agam, bgam, mu,Xij, b)
#
indrejtotal<-0 #initialize vector indrejtotal
probacetotal<-1 #initialize vector probacetotal
#
amostrasfin<-1000 # MCMC sample size after burn in and jumps
burnin<-1000 # burn in size
saltos<-10 # jumps size
AmostrasTotal<-burnin+amostrasfin*saltos # number of iterations to be run
set.seed(300)
#
library(compiler)
enableJIT(3)

#creo vettore di matrici Xplusj, j=1,...,J (una per paese)
library(MASS)
library(dplyr)
Xplus <- list()
for (j in 1:J){
  Xj = Xij[which(dataset$CNT==j),]
  Xplus[[j]] = solve((t(Xj)%*%Xj))%*%t(Xj)
}

#
for (int in 1:AmostrasTotal){
#
######### update Sj and Zij by a MH step of one j
#
J=J #numero totale paesi
j<-sample(1:J,1)
cat('\n', int, j, K)
candidato<-sample.tau.xi(dataset, Ij, Sj, Zij, tauj, K, Lk, theta.samp, j, alpha, beta, ni, mu, desvio,Xij,b)
paccept<-prob.accept(dataset, candidato[[8]], Sj, Zij, tauj, Lk, theta.samp, j, candidato[[3]], candidato[[2]], candidato[[6]], candidato[[7]], lambda, mu, agam, bgam, candidato[[4]], candidato[[5]], candidato[[10]], candidato[[11]],Xij,b)	[[1]]
probacetotal<-c(probacetotal,paccept)
aux2<-runif(1)
if (aux2<paccept){
	indrejtotal<-c(indrejtotal,0)
	Sj[j]<-candidato[[6]]
	tauj[(cumsum(Ij)[j]-Ij[j]+1):cumsum(Ij)[j]]<-rep(candidato[[6]],Ij[j])
	#modifichiamo noi le prossime 2 righe (prima non c'era l'if e il comando lo faceva eseguire sempre)
	if(sum(is.na(candidato[[7]]))==0)
	  Zij[(cumsum(Ij)[j]-Ij[j]+1):cumsum(Ij)[j]]<-candidato[[7]]
	while (length(table(Sj))<max(Sj)){ # exclude empty distributional clusters
		categr<-as.numeric(as.character(data.frame(table(Sj))[,1]))
		categd<-seq(1:length(table(Sj)))
		dif<-which(categr!=categd)
		for (i in 1:length(Sj)) if (Sj[i]>dif[1]) Sj[i]<-Sj[i]-1
		tauj<-numeric()
		for (n in 1:length(Sj)) tauj<-c(tauj,rep(Sj[n],Ij[n]))}	
	K<-length(table(Sj))
	Lk<-numeric(K)
	for (k in 1:K){
		observ<-which(tauj==k)
		Lk[k]<-length(table(Zij[observ]))
		while (Lk[k]<max(Zij[observ])){ # exclude empty observational clusters
			categr<-as.numeric(as.character(data.frame(table(Zij[observ]))[,1]))
			categd<-seq(1:Lk[k])
			dif<-which(categr!=categd)
			for (i in observ) if (Zij[i]>dif[1]) Zij[i]<-Zij[i]-1}}
	theta.samp<-posteriori.theta(dataset, Zij, tauj, lambda, agam, bgam, mu,Xij,b)
	Zij<-posteriori.Zij2(dataset, Zij, tauj, K, theta.samp, ni, mu, desvio, Lk, beta, j, Sj,Xij,b)
	Lk<-numeric(K)
	for (n in 1:K) Lk[n]<-length(table(Zij[which(tauj==n)]))
	theta.samp<-posteriori.theta(dataset, Zij, tauj, lambda, agam, bgam, mu,Xij,b)}
if (aux2>=paccept) indrejtotal<-c(indrejtotal,1)
#

bvec_totale <- matrix(nrow=J, ncol=4) 
for (j in 1:J) {
  muj=theta.samp[Zij[which(dataset$CNT==j)],1]
  sigmaj=theta.samp[Zij[which(dataset$CNT==j)],2]
  d <- sum(dataset$CNT == j)  # Conta quanti elementi hanno CNT == j
  inv_b=diag(4)*sigmab^(-2)
  inv_j=diag(sigmaj)
  
  
  # # Aggiunta di una regolarizzazione per evitare problemi di singolarità
  # lambda <- 1e-6  
  # 
  # # Matrice di precisione inversa regolarizzata
  # inv_b <- ginv(diag(4) * sigmab^(-2) + diag(lambda, 4))  
  # inv_j <- ginv(diag(sigmaj) + diag(lambda, d))  
  # 
  # # Calcola l'inversa regolarizzata usando SVD
  # XTX <- Xplus[[j]] %*% inv_j %*% t(Xplus[[j]]) + diag(lambda, nrow(Xplus[[j]]))  # Matrice regolarizzata
  # svd_XTX <- svd(XTX)
  # XTX_inv <- svd_XTX$v %*% diag(1 / (svd_XTX$d + lambda)) %*% t(svd_XTX$u)  # Pseudo-inversa regolarizzata
  # 
  # # Calcolo di bvec con l'inversa stabilizzata
  # mean_bvec <- XTX_inv %*% Xplus[[j]] %*% (dataset$PV1MATH[which(dataset$CNT==j)] - muj)
  # cov_bvec <- XTX_inv
  # 
 
  cov_bvec <- solve(inv_b + solve(Xplus[[j]]%*%inv_j%*%t(Xplus[[j]])))
  mean_bvec <- cov_bvec%*%solve(Xplus[[j]]%*%inv_j%*%t(Xplus[[j]]))%*%Xplus[[j]]%*%(dataset$PV1MATH[which(dataset$CNT==j)]-muj)
  bvec <- rmvnorm(1, mean_bvec, cov_bvec)  
  # Assegna i valori alla matrice finale
  bvec_totale[j, ] <- bvec
  
}
b <- bvec_totale
b1 <- b[,1]
b2 <- b[,2]
b3 <- b[,3]
b4 <- b[,4]
#



if (int>burnin & int%%saltos==0){
  write(Sj, file = paste(caminho, "Sj_simu2.txt", sep = ""))
  cat(Sj, "\n", file = paste(caminho, "Sj_simu1.txt", sep = ""), append = TRUE)
	cat('',theta.samp,file=paste(caminho,"theta_simu1.txt",sep=""),append=T)
	cat('',Zij,file=paste(caminho,"Zij_simu1.txt",sep=""),append=T)
	cat('',K,file=paste(caminho,"K_simu1.txt",sep=""),append=T)
	cat('',Lk,file=paste(caminho,"Lk_simu1.txt",sep=""),append=T)
	cat('',b1,file=paste(caminho,"beta1_simu1.txt",sep=""),append=T)
	cat('',b2,file=paste(caminho,"beta2_simu1.txt",sep=""),append=T)
	cat('',b3,file=paste(caminho,"beta3_simu1.txt",sep=""),append=T)
	cat('',b4,file=paste(caminho,"beta4_simu1.txt",sep=""),append=T)}
}
cat('',indrejtotal,file=paste(caminho,"indrejtotal_simu1.txt",sep=""),append=T)

unique(Sj)
unique(Zij)

### Riassegniamo le label dei cluster così che siano in ordine
cl <- rep(0,dim(dataset)[1])
for (i in 1:dim(dataset)[1]){
  cl[i] = Sj[dataset$CNT[i]]
}
cl <- as.factor(cl)
cl2 <- rep(0,dim(dataset)[1])
dataset$cl <- cl
for (i in 1:dim(dataset)[1]){
  if(cl[i]==1)
    cl2[i] = 8
  if(cl[i]==2)
    cl2[i] = 9
  if(cl[i]==3)
    cl2[i] = 4
  if(cl[i]==4)
    cl2[i] = 5
  if(cl[i]==5)
    cl2[i] = 3
  if(cl[i]==6)
    cl2[i] = 2
  if(cl[i]==7)
    cl2[i] = 1
  if(cl[i]==8)
    cl2[i] = 7
  if(cl[i]==9)
    cl2[i] = 6
}
unique(cl2)
cl <- as.factor(cl2)
data$cl <- cl
dataset$cl <- cl
 

### Distribution across countries
 library(dplyr)
library(ggplot2)
library(RColorBrewer)  
data <- data %>%
  mutate(cluster = cl)  # Convertiamo in fattore per colori discreti
distinct_colors <- c(
  "#00441B", "#1B7837", "#A6DBA0",  # 3 Verdi
  "#D0FF00", "#FFE119", "#F58231", # Giallo, Ocra, Arancione
  "#FF6347", "#D7191C", "#800020"   # Rosso Chiaro, Rosso, Bordeaux
)
ggplot(data, aes(x = as.factor(CNT), y = PV1MATH, color = factor(cluster))) +
  geom_jitter(width = 0.2, alpha = 0.7) +  
  scale_color_manual(values = distinct_colors[1:length(unique(data$cluster))]) +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Countries", y = "Students' average score", color = "Cluster") +
  ggtitle("Score distribution across countries")


### Mappa  
data2 <- read.csv("OECD_data.csv")
data2 = data2 %>% group_by(CNT) %>% filter(n()>3) %>% ungroup()
data2 = as.data.frame(data2)
clusters <- data.frame(CNT=data2$CNT,clust=dataset$cl)

library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)
library(readr)
library(rnaturalearthdata)

clusters <- clusters %>%
  mutate(CNT = case_when(
    CNT == "United States" ~ "United States of America",
    CNT == "B-S-J-Z (China)" ~ "China",
    CNT == "Chinese Taipei" ~ "Taiwan",
    CNT == "Slovak Republic" ~ "Slovakia",
    CNT == "Brunei Darussalam" ~ "Brunei",
    TRUE ~ CNT
  ))
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- world %>%
  left_join(clusters, by = c("name" = "CNT")) %>%
  filter(!sov_a3 %in% c("ATA", "GRL", "ATF"))  # Rimuove Antartide, Groenlandia e Terre Australi Francesi

n_clusters <- length(unique(na.omit(world$clust)))
distinct_colors <- c(
  "#00441B", "#1B7837", "#A6DBA0",  # 3 Verdi
  "#D0FF00", "#FFE119", "#F58231", # Giallo, Ocra, Arancione
  "#FF6347", "#D7191C", "#800020"   # Rosso Chiaro, Rosso, Bordeaux
)

ggplot() +
  geom_sf(data = world, aes(fill = factor(clust)), color = "whitesmoke") +
  scale_fill_manual(values = distinct_colors[1:n_clusters], name = "Cluster Index by Country", na.value = "grey80") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(title = "Cluster Map of Countries")



### Calcolo quanti sottocluster ci sono per ogni paese
tabella <- table(cl)
dfTemp <- data.frame(matrix(0, nrow = 38, ncol = 4))
for (i in 1:1223){
  cnt = dataset$CNT[i]
  cl = Zij[i]
  dfTemp[cnt,cl] <- dfTemp[cnt,cl] +1
}
dfTemp$Sj = Sj

### Convergenza MCMC

### Catena K -> numero cluster distribuzionali
K_chain <- scan(paste(caminho, "K_simu1.txt", sep=""))
burnin <- 2000  
K_chain <- K_chain[(burnin + 1):length(K_chain)]
plot(K_chain, type="l", col="blue", xlab="MCMC Iterations", 
     ylab="Number of Clusters", main="Number of clusters' Traceplot")

### Catena beta_1
beta1_chain <- scan(paste(caminho, "beta1_simu1.txt", sep=""))
num_paesi <- length(Sj)  
num_iterazioni <- length(beta1_chain) / num_paesi  
matrice_valori <- matrix(beta1_chain, nrow = num_iterazioni, ncol = num_paesi, byrow = TRUE) #ogni colonna un paese, ogni riga un'iterazione
valori_paese1 <- matrice_valori[,1]  
plot(valori_paese1, type="l", col="blue", xlab="Iterations MCMC", 
     ylab="beta1", main="Traceplot beta1 of Albania")

### Catena beta_2
beta2_chain <- scan(paste(caminho, "beta2_simu1.txt", sep=""))
num_paesi <- length(Sj)  
num_iterazioni <- length(beta2_chain) / num_paesi  
matrice_valori <- matrix(beta2_chain, nrow = num_iterazioni, ncol = num_paesi, byrow = TRUE) #ogni colonna un paese, ogni riga un'iterazione
valori_paese2 <- matrice_valori[,1]  
plot(valori_paese2, type="l", col="blue", xlab="Iterations MCMC", 
     ylab="beta2", main="Traceplot beta2 of Albania")

### Catena beta_3
beta3_chain <- scan(paste(caminho, "beta3_simu1.txt", sep=""))
num_paesi <- length(Sj)  
num_iterazioni <- length(beta3_chain) / num_paesi  
matrice_valori <- matrix(beta3_chain, nrow = num_iterazioni, ncol = num_paesi, byrow = TRUE) #ogni colonna un paese, ogni riga un'iterazione
valori_paese3 <- matrice_valori[,1]  
plot(valori_paese3, type="l", col="blue", xlab="Iterations MCMC", 
     ylab="beta3", main="Traceplot beta3 of Albania")

### Catena beta_4
beta4_chain <- scan(paste(caminho, "beta4_simu1.txt", sep=""))
num_paesi <- length(Sj)  
num_iterazioni <- length(beta4_chain) / num_paesi  
matrice_valori <- matrix(beta4_chain, nrow = num_iterazioni, ncol = num_paesi, byrow = TRUE) #ogni colonna un paese, ogni riga un'iterazione
valori_paese4 <- matrice_valori[,1]  
plot(valori_paese4, type="l", col="blue", xlab="Iterations MCMC", 
     ylab="beta4", main="Traceplot beta4 of Albania")

### Plot posterior-cluster
file_lines <- readLines(paste(caminho, "Sj_simu1.txt", sep = ""))
file_lines <- file_lines[-1]
writeLines(file_lines, paste(caminho, "Sj_simu1_corretto.txt", sep = ""))

Sj_matrix <- as.matrix(read.table(paste(caminho, "Sj_simu1_corretto.txt", sep = ""), header = FALSE))
num_iterazioni <- 999
coassignment_matrix <- diag(J)*num_iterazioni
for (t in 1:num_iterazioni) {
  for (i in 1:(J-1)) { 
    for (j in (i+1):J) { 
      if (Sj_matrix[t, i] == Sj_matrix[t, j]) {  
        coassignment_matrix[i, j] <- coassignment_matrix[i, j] + 1
        coassignment_matrix[j, i] <- coassignment_matrix[i, j]  
      }
    }
  }
}
coassignment_matrix <- coassignment_matrix / num_iterazioni
coassignment_matrix <- coassignment_matrix[nrow(coassignment_matrix):1, ]
J <- dim(coassignment_matrix)[2]
image(1:J, 1:J, coassignment_matrix, col = gray(seq(1, 0, length.out = 100)),
      xlab = "Paesi", ylab = "Paesi", main = "Matrice di Co-assegnazione", axes = FALSE)




### Caratterizzazione cluster 
library(fmsb)
covariate <- data[, !names(data) %in% c("cl","CNTSCHID","CNT","PRIVATESCH","sum_MATH1below","Y_BIN_MATH1","cluster","mean_ESCS_std","SCH_TESTED","Y_MATH1")]
cluster_means <- aggregate(covariate, by = list(dataset$cl), FUN = mean)
colnames(cluster_means)[1] <- "cl"
normalize <- function(x) (x - min(x)) / (max(x) - min(x) + 1e-9)  # Evita divisione per zero
cluster_means[, -1] <- as.data.frame(lapply(cluster_means[, -1], normalize))
num_clusters <- length(unique(dataset$cl))
par(mfrow = c(ceiling(num_clusters / 2), 2)) 
par(mfrow = c(1, 1))

i <- 1  # Cambia il numero per vedere un altro cluster
cluster_data <- cluster_means[i, -1] 
data_radar <- rbind(rep(1, ncol(cluster_data)), rep(0, ncol(cluster_data)), cluster_data)

radarchart(data_radar, 
           axistype = 2,
           pcol = "blue", pfcol = rgb(0.2, 0.5, 0.5, 0.3), 
           title = paste("Radarchart cluster 1", cluster_means$Cluster[i]))

i <- 9  # Cambia il numero per vedere un altro cluster
cluster_data <- cluster_means[i, -1] 
data_radar <- rbind(rep(1, ncol(cluster_data)), rep(0, ncol(cluster_data)), cluster_data)

radarchart(data_radar, 
           axistype = 2,
           pcol = "blue", pfcol = rgb(0.2, 0.5, 0.5, 0.3), 
           title = paste("Radarchart cluster 9", cluster_means$Cluster[i]))


i <- 4  # Cambia il numero per vedere un altro cluster
cluster_data <- cluster_means[i, -1] 
data_radar <- rbind(rep(1, ncol(cluster_data)), rep(0, ncol(cluster_data)), cluster_data)

radarchart(data_radar, 
           axistype = 2,
           pcol = "blue", pfcol = rgb(0.2, 0.5, 0.5, 0.3), 
           title = paste("Radarchart cluster 4", cluster_means$Cluster[i]))

### Credible intervals

## Deduco nomi paesi
country_names <- data2$CNT
country_numbers <- dataset$CNT
df_nomi <- data.frame(Number = country_numbers, Country = country_names)
df_nomi <- unique(df_nomi)
df_nomi <- df_nomi[order(df_nomi$Number), ]

## Costruisco credible intervals dei BETA2 
beta_samples2 <- matrix(beta2_chain, nrow = num_iterazioni, ncol = num_paesi, byrow = TRUE)
credible_intervals <- apply(beta_samples2, 2, credible_interval, level = 0.95)
colnames(beta_samples2) <- df_nomi$Country
colnames(credible_intervals) <- df_nomi$Country 

# Convertiamo in un dataframe per ggplot
credible_df <- data.frame(
  State = colnames(beta_samples2),
  Lower = credible_intervals[1, ],
  Upper = credible_intervals[2, ],
  Mean = colMeans(beta_samples2)
)
country_clusters <- dataset %>%
  group_by(`CNT`) %>%  # Raggruppiamo per numero di country
  summarise(cluster = unique(cl)) %>%  # Prendiamo il valore unico di cluster per ogni paese
  arrange(`CNT`) %>%  # Ordiniamo per numero di country
  pull(cluster)  # Convertiamo in vettore
credible_df$Cluster <- country_clusters

ggplot(credible_df, aes(x = reorder(State, Mean), y = Mean, color = factor(Cluster))) +
  geom_point(size = 3) +  # Punto per il valore medio
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +  # Intervalli credibili
  theme_minimal() +
  scale_color_manual(values = distinct_colors) +  # Assegniamo i colori ai cluster
  labs(x = "State", y = "Beta2 Estimate", title = "Credible Intervals for Beta2 by State", color = "Cluster") +
  coord_flip()  # Migliora leggibilità

## Costruisco credible intervals dei BETA1 
credible_interval <- function(samples, level = 0.95) {
  lower <- quantile(samples, probs = (1 - level) / 2)
  upper <- quantile(samples, probs = 1 - (1 - level) / 2)
  return(c(lower, upper))
}

beta_samples1 <- matrix(beta1_chain, nrow = num_iterazioni, ncol = num_paesi, byrow = TRUE)
credible_intervals <- apply(beta_samples1, 2, credible_interval, level = 0.95)
colnames(beta_samples1) <- df_nomi$Country
colnames(credible_intervals) <- df_nomi$Country 

# Convertiamo in un dataframe per ggplot
credible_df <- data.frame(
  State = colnames(beta_samples1),
  Lower = credible_intervals[1, ],
  Upper = credible_intervals[2, ],
  Mean = colMeans(beta_samples1)
)
country_clusters <- dataset %>%
  group_by(`CNT`) %>%  # Raggruppiamo per numero di country
  summarise(cluster = unique(cl)) %>%  # Prendiamo il valore unico di cluster per ogni paese
  arrange(`CNT`) %>%  # Ordiniamo per numero di country
  pull(cluster)  # Convertiamo in vettore
credible_df$Cluster <- country_clusters

ggplot(credible_df, aes(x = reorder(State, Mean), y = Mean, color = factor(Cluster))) +
  geom_point(size = 3) +  # Punto per il valore medio
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +  # Intervalli credibili
  theme_minimal() +
  scale_color_manual(values = distinct_colors) +  # Assegniamo i colori ai cluster
  labs(x = "State", y = "Beta 1 Estimate", title = "Credible Intervals for Beta 1 by State", color = "Cluster") +
  coord_flip()  # Migliora leggibilità

## Costruisco credible intervals dei BETA3 
beta_samples3 <- matrix(beta3_chain, nrow = num_iterazioni, ncol = num_paesi, byrow = TRUE)
credible_intervals <- apply(beta_samples3, 2, credible_interval, level = 0.95)
colnames(beta_samples3) <- df_nomi$Country
colnames(credible_intervals) <- df_nomi$Country 

# Convertiamo in un dataframe per ggplot
credible_df <- data.frame(
  State = colnames(beta_samples3),
  Lower = credible_intervals[1, ],
  Upper = credible_intervals[2, ],
  Mean = colMeans(beta_samples3)
)
country_clusters <- dataset %>%
  group_by(`CNT`) %>%  # Raggruppiamo per numero di country
  summarise(cluster = unique(cl)) %>%  # Prendiamo il valore unico di cluster per ogni paese
  arrange(`CNT`) %>%  # Ordiniamo per numero di country
  pull(cluster)  # Convertiamo in vettore
credible_df$Cluster <- country_clusters

ggplot(credible_df, aes(x = reorder(State, Mean), y = Mean, color = factor(Cluster))) +
  geom_point(size = 3) +  # Punto per il valore medio
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +  # Intervalli credibili
  theme_minimal() +
  scale_color_manual(values = distinct_colors) +  # Assegniamo i colori ai cluster
  labs(x = "State", y = "Beta 3 Estimate", title = "Credible Intervals for Beta 3 by State", color = "Cluster") +
  coord_flip()  # Migliora leggibilità

## Costruisco credible intervals dei BETA4 
beta_samples4 <- matrix(beta4_chain, nrow = num_iterazioni, ncol = num_paesi, byrow = TRUE)
credible_intervals <- apply(beta_samples4, 2, credible_interval, level = 0.95)
colnames(beta_samples4) <- df_nomi$Country
colnames(credible_intervals) <- df_nomi$Country 

# Convertiamo in un dataframe per ggplot
credible_df <- data.frame(
  State = colnames(beta_samples4),
  Lower = credible_intervals[1, ],
  Upper = credible_intervals[2, ],
  Mean = colMeans(beta_samples4)
)
country_clusters <- dataset %>%
  group_by(`CNT`) %>%  # Raggruppiamo per numero di country
  summarise(cluster = unique(cl)) %>%  # Prendiamo il valore unico di cluster per ogni paese
  arrange(`CNT`) %>%  # Ordiniamo per numero di country
  pull(cluster)  # Convertiamo in vettore
credible_df$Cluster <- country_clusters

ggplot(credible_df, aes(x = reorder(State, Mean), y = Mean, color = factor(Cluster))) +
  geom_point(size = 3) +  # Punto per il valore medio
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +  # Intervalli credibili
  theme_minimal() +
  scale_color_manual(values = distinct_colors) +  # Assegniamo i colori ai cluster
  labs(x = "State", y = "Beta 4 Estimate", title = "Credible Intervals for Beta 4 by State", color = "Cluster") +
  coord_flip()  # Migliora leggibilità


