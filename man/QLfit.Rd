% File src/library/QLSpline/man/QLfit.Rd
     \name{QLfit}
     \alias{QLfit}
     \title{Fit the quasi-likelihood model to matrix of RNA-seq expression count data}
     \description{
       Analyze RNA-seq expression count data using the methods detailed in Lund and Nettleton (2011).  This is the main function of the QLSpline package.}
     \usage{
    QLfit(dat,design.list,test.mat=NULL,log.offset=NULL,spline.df=NULL,Plot=FALSE)
     }
     \arguments{
      \item{dat}{RNA-seq data matrix of integer expression counts.  Each row contains observations from a single gene. Each column contains observations from a single experimental unit.}
      \item{design.list}{List of design matrices for the full model and reduced model(s). The first element of \code{design.list} must describe the overall full model, as this design is used to compute deviance residuals for estimating dispersion.  One factor designs may be specified as vectors. The number of rows in each design matrix (or the length of each design vector) must be \code{ncol(dat)}.}
	\item{test.mat}{T by 2 matrix dictating which designs are to be compared, where T is the total number of desired hypothesis tests for each gene.  Each row contains two integers, which provide the indices within \code{design.list} of the two designs to be compared.  If \code{test.mat} is not specified, the default is compare the first design (the full model) to each of the other designs provided in \code{design.list}.}
	\item{log.offset}{A vector of log-scale, additive factors used to adjust estimated log-scale means for differences in library sizes across samples.  Commonly used offsets include,\code{log.offset=log(apply(dat,2,sum))} or \code{log.offset=log(apply(dat,2,quantile,.75))}.  The default setting makes no adjustment for library sizes (i.e. log.offset=0).}
	\item{spline.df}{Optional.  User may specify the degrees of freedom to use when fitting a cubic spline to log-scale(estimated dispersion) versus the log(average count).  Default uses cross-validation in \code{sreg} function to pick appropriate degrees of freedom.} 	
	\item{Plot}{logical.  If TRUE, the estimated dispersion versus the average count are plotted on a log-scale with the corresponding cubic spline fit overlaid. }
     }

\value{list containing:
	\item{"LRT"}{matrix providing unadjusted likelihood ratio test statistics.  Each column contains statistics from a single hypothesis test, applied separately to each gene.}
	\item{"phi.hat.dev"}{vector providing unshrunken, deviance-based estimates of dispersion (phi) for each gene.}
%	\item{"num.df"}{vector of numerator degrees of freedom for each test, computed as the difference in the number of fitted parameters between the full and reduced models for each test.}
%	\item{"den.df"}{denominator degrees of freedom. Equal to the number of samples minus the number of fitted parameters in the model specified by the first element of \code{design.list}}
	\item{"P.values"}{list of matrices providing p-values.  The i^th element of \code{pvals} corresponds to the hypothesis test assigned in the i^th row of \code{test.mat}.
If there is only one test for each gene, \code{pvals} is returned as a matrix, rather than a list.  Each matrix contains six columns, one for each of the following approaches: QL.F, QL.Chi, QLShrink.F, QLShrink.Chi, QLSpline.F, and QLSpline.Chi.  Note that column names for QLShrink.F and QLSpline.F methods contain number of additional denominator d.f. gained from shrinking dispersion estimates.}
\item{"Q.values"}{list of matrices providing q-values.  The i^th element of \code{qvals} corresponds to the hypothesis test assigned in the i^th row of \code{test.mat}.
If there is only one test for each gene, \code{qvals} is returned as a matrix, rather than a list.  Each matrix contains six columns, one for each of the following approaches: QL.F, QL.Chi, QLShrink.F, QLShrink.Chi, QLSpline.F, and QLSpline.Chi.  Q-values are computed using the methods of Nettleton et al. (2006) JABES 11, 337-356.}
\item{"m0"}{matrix providing estimated number of true null hypotheses for each test(row) under each of the six methods(columns). \code{m0} values are computed using the methods of Nettleton et al. (2006) JABES 11, 337-356.}
}

\author{Steve Lund \email{lunds@iastate.edu}}



\examples{  
### Create example data set from a gammma-Poisson model.  
trt<-rep(1:2,each=5);n<-length(trt); n.genes<-2000

#### Simulate means for 2000 genes
sim.mn<-matrix(exp(rnorm(n.genes,2.5,2)),n.genes,2)

#### make first 1000 genes DE
sim.mn[1:(.5*n.genes),1]<-exp(rnorm(.5*n.genes,2.5,2))

### Simulate dispersion parameters
phi<-rgamma(n.genes,2,.7)+1

b<-1/(phi-1); a<-sim.mn[,trt]*b

### Simulate library size factors
offset<-2^(rnorm(n,0,.5))
a<-t(t(a)*offset)
lambda<-matrix(rgamma(n.genes*n,a,rep(b,n)),n.genes,n)
simdat<-matrix(rpois(n.genes*n,lambda),n.genes,n)

### Keep genes with at least 10 total counts
simdat<-simdat[rowSums(simdat)>9,]

### Create list of designs describing model under null and alternative hypotheses
design.list<-vector("list",2)
design.list[[1]]<-model.matrix(~as.factor(trt))  #This also could have just been ``trt''.
design.list[[2]]<-rep(1,length(trt))

log.offset<-log(apply(simdat,2,quantile,.75))

### Analyze using QL, QLShrink, and QLSpline methods
results<-QLfit(simdat, design.list,log.offset=log.offset)

### How many significant genes at FDR=.05?
apply(results$Q.values<.05,2,sum)

### Indexes for Top 100 most significant genes from QLShrink method
order(results$P.values[,5])[1:100]

}
     \keyword{RNA-seq,quasi-likelihood, differential expression}



