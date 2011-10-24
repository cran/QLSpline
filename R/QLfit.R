
QLfit<-function(dat,design.list,test.mat=NULL,log.offset=NULL,spline.df=NULL,Plot=FALSE){

### First element of design.list should pertain to overall full model.  This is the design used to obtain phi.hat.dev.

### Evaluate the deviance under each design provided in design.list
deviance.list<-vector("list",length(design.list))

p<-NULL; n<-ncol(dat)
for(jj in 1:length(design.list)){
	design<-design.list[[jj]]

####  collapse a one factor design matrix into vector for increased speed
if(is.matrix(design)) {
if(ncol(design)==1) design<-as.vector(design)
if(prod(sort(unique(as.vector(design)))==c(0,1))==1){   ### Is design matrix made up of only 1's and 0's?
if(prod(design[,1]==1)==1){ 	 				 ### First column all 1s for intercept?                
if(max(rowSums(design))<3){ design<-as.vector(design%*%c(1,1:(ncol(design)-1)))   ### Does each row include at most one additional 1? Then it's a one factor design
} } } }


	if(is.matrix(design)) {
	### For multi-factor designs, the first column of each element (matrix) of design.list should be a column of 1's, pertaining to the intercept. 
		##offset should be given on log scale here
		apply.glm<-function(y,design,log.offset=NULL) glm(y~design[,-1],family="poisson",offset=log.offset)$deviance 
		p<-c(p,ncol(design))

### Screen for basic errors
		if(prod(design[,1]==1)!=1) stop(paste("The first column of matrix in element",jj,"of 'design.list' is not a column of 1s for the intercept. Please include intercept."))
		if(nrow(design)!=n) stop(paste("Element",jj,"in 'design.list' has",nrow(design),"rows. 
Design matrices must have",n,"rows (to match number of columns in data)."))
		if(p[jj]>p[1]) stop(paste("Full model design must be element in 'design.list'. 
'p' for element",jj,"is larger than 'p' for first element, 
indicating first element does not correspond to full model design."))
		deviance.list[[jj]]<-apply(dat,1,apply.glm,design=design,log.offset=log.offset)
	}

### For one factor designs, MLE's for Poisson GLM have simple closed form, so we can avoid one-gene-at-a-time GLM fitting.
	if(is.vector(design)) {
		p<-c(p,length(unique(design)))
### Screen for basic errors
		if(p[jj]>p[1]) stop(paste("Full model design must be element in 'design.list'. 
'p' for element",jj,"is larger than 'p' for first element, 
indicating first element does not correspond to full model design."))
		if(length(design)!=n) stop(paste("Element",jj,"in 'design.list' has length",length(design),". 
Design vectors must have length",n,"(to match number of columns in data)."))

		##offset should NOT be given on log scale here
		offset<-rep(1,length(design)); if(!is.null(log.offset)) offset<-exp(log.offset)

		dat<-as.matrix(dat)
		means<-dat

		for(i in 1:p[jj]) means[,design==unique(design)[i]]<-rowSums(dat[,design==unique(design)[i]])/sum(offset[design==unique(design)[i]])
		means<-t(t(means)*offset)

		### 0's require special attention since 0^0=1, but 0*log(0)=NaN
		deviance<-means-dat
		deviance[dat!=0]<-deviance[dat!=0]+(dat*log(dat/means))[dat!=0]
		deviance.list[[jj]]<-2*rowSums(deviance)
	}
}

### Add ons:

###  confirm that that first design matrix is largest (i.e. corresponds to full model)
###  confirm that design matrices include intercept column


### If not otherwise specified, compare each design matrix to the first design matrix, 
### which should be the full design matrix
if(is.null(test.mat)){ test.mat<-cbind(1,2:length(design.list))
	rownames(test.mat)<-paste(1,2:length(design.list))  }

LRT<-NULL;num.df<-NULL
for(i in 1:nrow(test.mat)){
	i1<-test.mat[i,1]; i2<-test.mat[i,2]
	num.df<-c(num.df,abs(p[i2]-p[i1]))
	LRT<-cbind(LRT,abs((deviance.list[[i2]]-deviance.list[[i1]])/(p[i2]-p[i1])))
}
colnames(LRT)<-rownames(test.mat)

den.df<-(n-p[1])

phi.hat.dev<-deviance.list[[1]]/den.df


### We use Smyth's (2004) approach from LIMMA to estimate parameters for prior distributions 
### of gene specific dispersion estimates. The function below also provides resulting 
### shrunken point estimates of dispersion 

#### Code for implementing Smyth's approach begins here ####
	shrink.phi<-function(phi.hat,den.df){
	phi.hat[phi.hat<=0]<-min(phi.hat[phi.hat>0])
	z<-log(phi.hat); z[z==Inf]<-max(z[z!=Inf]); z[z==-Inf]<-min(z[z!=-Inf]);mnz<-mean(z)

	## solve for d0 and phi0
	d0arg<-var(z)-trigamma((den.df)/2)
	if(d0arg>0){
		dif<-function(x,y) abs(trigamma(x)-y)
		inverse.trigamma<-function(y) optimize(dif,interval=c(0,10000),y=y)$minimum
		d0<-2*inverse.trigamma(d0arg)
		phi0<-exp(mnz-digamma((den.df)/2)+digamma(d0/2)- log(d0/(den.df)))

		## compute shrunken phi's
		phi.shrink<-((den.df)*phi.hat+d0*phi0)/(den.df+d0-2)  }
	else{phi.shrink<-rep(exp(mnz),length(z)); d0<-Inf; phi0<-exp(mnz) }
	return(list(phi.shrink=phi.shrink,d0=d0,phi0=phi0))  }
#### Code for implementing Smyth's approach ends here ####


pval<-vector("list",ncol(LRT)); names(pval)<-colnames(LRT)

phi.hat.dev[phi.hat.dev<0]<-min(phi.hat.dev[phi.hat.dev>0])
phi.hat.dev2<-phi.hat.dev; phi.hat.dev2[phi.hat.dev<1]<-1

shrink<-shrink.phi(phi.hat.dev,den.df)
phi.shrink<-shrink[[1]]; est.d0<-shrink[[2]]; est.s02<-shrink[[3]]
if(est.d0==Inf) est.d0<-10

phi.shrink[phi.shrink<1]<-1

### Fit cubic spline to capture relationship between log(y.bar) and log(phi.hat.dev)
rowsums.dat<-rowSums(dat[,])
y<-log(phi.hat.dev); y[y==-Inf]<-min(y[y!=-Inf]); y[y==Inf]<-max(y[y!=Inf])
if(is.null(spline.df)) spline.fit<-sreg(log(rowsums.dat/ncol(dat)), y)
else spline.fit<-sreg(log(rowsums.dat), y,df=spline.df)

### If desired, plot resulting cubic spline fit
if(Plot){
	windows(); par(mai=c(1,1.2,1,.2))
	plot(log(rowsums.dat/ncol(dat)),y,xlab=expression(log(bar(Y)[phantom()%.%phantom()]*phantom()[phantom()%.%phantom()])),
ylab=expression(log(hat(Phi))),main="Estimated Dispersion
 versus Average Count",pch=1,cex.lab=2,cex.axis=2,cex.main=2)
	lines(sort(log(rowsums.dat/ncol(dat))),spline.fit$fitted.values[order(rowsums.dat)],col=2,lwd=3)
	legend("bottomright",legend=paste("Fitted Spline with ",signif(spline.fit$eff.df,2),"df"),lwd=3,col=2,cex=1.5)
}

### Obtain estimate for prior degrees of freedom after adjusting for cubic spline fit
y2<-phi.hat.dev/exp(spline.fit$fitted.values)
shrink<-shrink.phi(y2,den.df)
D0.new<-shrink[[2]]; 
phi0.new<-shrink[[3]]

phi.spline.new<-(D0.new*exp(spline.fit$fitted.values)/phi0.new+(den.df)*phi.hat.dev)/(D0.new+den.df-2)
phi.spline.new[phi.spline.new<1]<-1

#### Compute p-values for each hypothesis test using each of the six approaches.
for(i in 1:ncol(LRT)){

pval[[i]]<-cbind(1-pf(LRT[,i]/phi.hat.dev2,num.df[i],den.df), 1-pf(LRT[,i]/phi.hat.dev2,num.df[i],Inf),
1-pf(LRT[,i]/phi.shrink,num.df[i],est.d0+den.df), 1-pf(LRT[,i]/phi.shrink,num.df[i],Inf), 
1-pf(LRT[,i]/phi.spline.new,num.df[i],D0.new+den.df),1-pf(LRT[,i]/phi.spline.new,num.df[i],Inf))

colnames(pval[[i]])<-c("QL.F","QL.Chi",paste("QLShrink.F d0=",round(est.d0,1),sep=""),"QLShrink.Chi",paste("QLSpline.F d0=",round(D0.new,1),sep=""),"QLSpline.Chi") }

estimate.m0<-function(p, B = 20){

#This function estimates the number of true null hypotheses given a vector of p-values
#using the method of Nettleton et al. (2006) JABES 11, 337-356.
#The estimate obtained is identical to the estimate obtained by the iterative
#procedure described by Mosig et al. Genetics 157:1683-1698.
#The number of p-values falling into B equally sized bins are counted.
#The count of each bin is compared to the average of all the bin counts associated
#with the current bins and all bins to its right.  Working from left to right, 
#the first bin that has a count less than or equal to the average is identified.
#That average is multiplied by the total number of bins to obtain an estimate of m0, 
#the number of tests for which the null hypothesis is true.
#
#Author: Dan Nettleton
  m <- length(p)
  m0 <- m
  bin <- c(-0.1, (1:B)/B)
  bin.counts=rep(0,B)
  for(i in 1:B){
    bin.counts[i]=sum((p>bin[i])&(p<=bin[i+1]))
  }
  tail.means <- rev(cumsum(rev(bin.counts))/(1:B))
  temp <- bin.counts - tail.means
  index <- min((1:B)[temp <= 0])
  m0 <- B * tail.means[index]
  return(m0)
}

jabes.q<-function(p,B=20){
#
#This function computes q-values using the approach of Nettleton et al.
#(2006) JABES 11, 337-356.
#
#Author: Dan Nettleton
#
  
  m = length(p)
  m0=estimate.m0(p,B)
  k = 1:m
  ord = order(p)
  p[ord] = (p[ord] * m0)/(1:m)
  qval = p
  qval[ord]=rev(cummin(rev(qval[ord])))
  return(qval)
}

#### Use p-values to obtain corresponding q-values 
#### and estimated proportions of null genes (m0)
qval<-pval; m0<-NULL
for(ii in 1:length(qval)){
	M0<-NULL; Qval<-qval[[ii]]
	for(jj in 1:ncol(Qval)){
		M0<-c(M0,estimate.m0(Qval[!is.na(Qval[,jj]),jj]))
		qval[[ii]][!is.na(Qval[,jj]),jj]=jabes.q(Qval[!is.na(Qval[,jj]),jj])
	}
	m0<-rbind(m0,M0)
}

colnames(m0)<-colnames(qval[[1]]); rownames(m0)<-names(pval)

### If there's only one test, use matrices instead of lists
if(length(pval)==1){pval<-pval[[1]];qval<-qval[[1]];m0<-m0[1,];names(m0)<-colnames(qval)}



#return(list(LRT=LRT,phi.hat.dev=phi.hat.dev,num.df=num.df,den.df=(n-p[1])))


return(list(LRT=LRT,phi.hat.dev=phi.hat.dev,P.values=pval,Q.values=qval,m0=m0))
}
