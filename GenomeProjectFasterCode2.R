#######################################################################################
## load the functions
#source('~/REU12/src/R/commonPre-OptimizationLibrary.R')
source("D:/cygwin/home/Administrator/REU12/src/R/commonPre-OptimizationLibrary.R") #work on Claire's Laptop
#source('../Desktop/REU12/src/R/estimateParaWithPhi.R')
source('D:/cygwin/home/Administrator/REU12/src/R/estimateParaWithPhi.R')
###source('~/REU12/src/R/scratch.R')
source('D:/cygwin/home/Administrator/REU12/scratch/chaij/loadlibrary.R')
#source('C:/cygwin/home/owner/REU12/scratch/liwe01/loadlibrary.R')
#source('C:/cygwin/home/owner/REU12/scratch/chaij/multinomial.R')
source('D:/cygwin/home/Administrator/REU12/scratch/chaij/multinomial.R')
#load the simulated codon counts variable saved in R workspace
load('D:/cygwin/home/Administrator/REU12/data/input_files/codon_counts/S.cerevisiae.S288c.simulated_1.RData')
######################################################################################################
elongtimeFile <- "D:/cygwin/home/Administrator/REU12/data/input_files/elong_times/S.cerevisiae/S.cere.elong.times.beyer.whole.genome.tsv"
mutationFile <- "D:/cygwin/home/Administrator/REU12/data/input_files/mutation_rates/S.cerevisiae/S.cere.mut.rates.beyer.whole.genome.tsv"
fastaFile <- "D:/cygwin/home/Administrator/REU12/data/input_files/FASTA/S.cerevisiae/S.cerevisiae.S288c.fasta"

#create the cCountsObj
cCountsObj <- createCodonCount(fastaFile)
elongTimeObj <- read_elong_time_file_and_create_aaList_and_nCodons(elongtimeFile)
elongTimeTable<-elongTimeObj$elongTime#elongation times for each synonymous codon of each aa
aaList <- elongTimeObj$aaList#list of all the aa's
nAA <- length(aaList)#length of the aaList(19)
nCodons <-elongTimeObj$nCodons#number of synonmyous codons per each aa
######################################################################################################
#This is a function that reads .tsv file and returns a vector of elongation time /mutation rates for each synonymous codon of each aa (maxAA=19)
get_ElongTime_MutRate_List <- function(filename)
{
  elongTimeObj <- read_elong_time_file_and_create_aaList_and_nCodons(filename) #elongTimeObj is a list with 3 components: a table consisting of elongation times, aalist, and nCodons 
  elongTimeTable<-elongTimeObj$elongTime #extract the elongTime table from the list above, first element of Obj
  nCodons <- elongTimeObj[[3]] #List of amino acids with their coding codon numbers in the file
  aaList <- elongTimeObj$aaList # list of amino acids in the file
  nAA <- length(aaList) #number of different amino acids(19)
  delta.t <- vector("list",length=nAA) # elongation time differences to return in vector form
  names(delta.t) <- aaList #change names of delta.t from numbers to corresponding aa's
  for(aa in aaList)
    delta.t[[aa]] <- elongTimeTable[elongTimeTable[,"aa"]==aa,"value"] #for each amino acid, extract times
  return(delta.t)
}

UnscaledDeltaTs <- get_ElongTime_MutRate_List(elongtimeFile)
DeltaT <- function(UnscaledDeltaTs){                     
  elongDiffObj <- lapply(UnscaledDeltaTs, function(x){x-min(x)})
  return(elongDiffObj)
}

mutRate <- get_ElongTime_MutRate_List(mutationFile)
delta.t <- DeltaT(UnscaledDeltaTs)
################################################################################


#ORF-names of ORF extracted from cCountsObj
ORF <-cCountsObj[,"ORF"]
ORFnamesListLength <-length(unique(ORF)) #number of differnt genes
ORFnamesUnique <- unique(ORF) #all different gene names without repeats
# ORFList <- vector("list", length = ORFnamesListLength)
#names(ORFList) <- unique(ORF)
reorderedcCountsObj <- cCountsObj[with(cCountsObj, order(ORF,aa)),]

simu.list <-simulated.codoncounts.list
simu.cCountsObj <- ldply(simu.list) ##
simu.ORF<-simu.cCountsObj[,"ORF"]
simu.ORFnamesListLength <-length(unique(simu.ORF))
simu.ORFnamesUnique <- unique(simu.ORF)
colnames(simu.cCountsObj)[1]<-"aa"
reorderedsimu.cCountsObj <- simu.cCountsObj[with(simu.cCountsObj, order(ORF,aa)),]


extractSameORF <- function(anORF){
  #out.list <- c(out.list,reorderedcCountsObj[reorderedcCountsObj[,"ORF"]==anORF,-2])
  theSameORF <- as.matrix.cast_matrix(reorderedcCountsObj[reorderedcCountsObj[,"ORF"]==anORF,-2]) 
  #names(out.list) <- anORF
  #tmp <- as.data.frame.cast_matrix(reorderedcCountsObj[reorderedcCountsObj[,"ORF"]==anORF,])
}
theSameORFList <- lapply(ORFnamesUnique,extractSameORF)
#tmp <- dataobj[[1]]
#tmp <- as.matrix.cast_matrix(tmp)
#rownames(tmp) <- NULL
names(theSameORFList) <- ORFnamesUnique
simu.theSameORFList <-lapply(simu.ORFnamesUnique,extractSameORF)
names(simu.theSameORFList)<-simu.ORFnamesUnique


getNumeratorAA<- function(mutRateOfaa,delta.tOfaa,z=1)
{
  return(mutRateOfaa*exp(-z*delta.tOfaa))
}


#This is a function that attains the likelihood for a particular gene; ln if/else statement included in this version 
#definition of terms and comments are same as in function above
getLikelihoodForAnORF <- function(theSameORF,z=1){ 

  
  theFinalResult <- 1  
  for (aa in aaList)
  {
    delta.tOfaa <- delta.t[[aa]] 
    mutRateOfaa<- mutRate[[aa]]
    nC <- nCodons[[aa]]
    #cCountsOfaa<-theSameORF[theSameORF[,"aa"] == aa,][2:(2+nC-1)] #number of observations of the synonymous codons in a particular amino acid
    cCountsOfaa<-theSameORF[theSameORF[,"aa"] == aa,][2:(2+nC-1)]
    #print(cCountsOfaa)
    NumList <- getNumeratorAA(mutRateOfaa,delta.tOfaa,z)
    NumListSum <-sum(NumList)
    total.cCountsPerAA <- sum(cCountsOfaa)#total number of observed synonymous codon counts for a particular aa
    NumListRaisedto.cCounts <- as.numeric(NumList^cCountsOfaa)
    NumListRaisedto.cCountsProd<-prod(NumListRaisedto.cCounts)
    result <- NumListRaisedto.cCountsProd * ((1/NumListSum)^total.cCountsPerAA)
    theFinalResult <- theFinalResult * result    
  }
  return (theFinalResult)
}



## Log version of the function above
getLnLikelihoodForAnORF <- function(theSameORF,z=1){ 
  
  
  theFinalResult <- 1  
  for (aa in aaList)
  {
    delta.tOfaa <- delta.t[[aa]] 
    mutRateOfaa<- mutRate[[aa]]
    nC <- nCodons[[aa]]
    cCountsOfaa<-theSameORF[theSameORF[,"aa"] == aa,][2:(2+nC-1)] #number of observations of the synonymous codons in a particular amino acid
    ##NumList <- getNumeratorAA(mutRateOfaa,delta.tOfaa,z)
    LnNumList <- log(mutRateOfaa)-delta.tOfaa*z
    NumListSum <-sum(exp(LnNumList))
    total.cCountsPerAA <- sum(cCountsOfaa)#total number of observed synonymous codon counts for a particular aa
    lnNumListRaisedto.cCounts<- as.numeric(cCountsOfaa*LnNumList)
    lnNumListRaisedto.cCountsProd<-sum(lnNumListRaisedto.cCounts)  
    result <- lnNumListRaisedto.cCountsProd - total.cCountsPerAA * log(NumListSum)
    theFinalResult <- theFinalResult + result    
  }
 
  return (theFinalResult)
  
}

##########################################################################
#This function is a modified version of the multinomial function built to attain the multinomial coefficient for a particular gene 
getMultinomForAnORF <- function(theSameORF){
  theFinalResult <- 1  
  for (aa in aaList)
  {
    nC <- nCodons[[aa]]#the number of types of synonymous codons for a particular aa
    cCountsOfaa<-theSameORF[theSameORF[,"aa"] == aa,][2:(2+nC-1)]#number of observations of the synonymous codons in a particular amino acid
    MultinomPerAA <- multinom(cCountsOfaa)
#     totalOf.cCountsPerAA <- sum(cCountsOfaa)#total number of observed synonymous codon counts for a particular aa
#     NumMultinom <- factorial(totalOf.cCountsPerAA)#numerator of the multinomial coefficient equation
#     DenomMultinomPreProduct <- factorial(cCountsOfaa)#denominator of the multinomial coefficient prior to multiplying across all synonymous codons of an aa
#     DenomMultinom <- prod(DenomMultinomPreProduct)#denominator of the multinomial coefficient after multiplying across synonmyous codons of an aa
#     MultinomPreProduct <- NumMultinom/DenomMultinom#multinomial coefficient for each aa prior to multiplying across all aa's     
    theFinalResult <- MultinomPerAA * theFinalResult#multinomial coefficient for particular gene sequence including all aa's
  }
  return(theFinalResult)#produces the multinomial coefficient for a specific gene
}

#This function is a modified version of the multinomial function built to attain the log of the multinomial coefficient for a particular gene rather than producing a list of non-log coefficients for all genes
getLnMultinomForAnORF <- function(theSameORF){ 
  theFinalResult <- 1
  for (aa in aaList)
  {
    nC <- nCodons[[aa]]#the number of types of synonymous codons for a particular aa
    cCountsOfaa<-theSameORF[theSameORF[,"aa"] == aa,][2:(2+nC-1)]#number of observations of the synonymous codons in a particular amino acid
    MultinomPerAA <- lmultinom(cCountsOfaa)
    #     totalOf.cCountsPerAA <- sum(cCountsOfaa)#total number of observed synonymous codon counts for a particular aa
    #     NumMultinom <- factorial(totalOf.cCountsPerAA)#numerator of the multinomial coefficient equation
    #     DenomMultinomPreProduct <- factorial(cCountsOfaa)#denominator of the multinomial coefficient prior to multiplying across all synonymous codons of an aa
    #     DenomMultinom <- prod(DenomMultinomPreProduct)#denominator of the multinomial coefficient after multiplying across synonmyous codons of an aa
    #     MultinomPreProduct <- NumMultinom/DenomMultinom#multinomial coefficient for each aa prior to multiplying across all aa's     
    theFinalResult <- MultinomPerAA + theFinalResult#multinomial coefficient for particular gene sequence including all aa's
  }
  return(theFinalResult)#produces the log of the multinomial coefficient for a specific gene
}
##########################################################################
#We now want to take logs of both the multinomial coefficient and likelihood for a particular gene and add them(thus attaining the product), and finally exponentiate them to attain the non-log value of said product
#This section is a function adding the logs of the multinomial coefficient and likelihood for a specific gene together
logLikelihoodMultinomProduct <- function(theSameORF, z){
  getLnMultinomForAnORF(theSameORF) + getLnLikelihoodForAnORF(theSameORF,z)
  }

logUnscaledPosterior <- function(theSameORF,z,meanlog,sdlog,log = TRUE){
  logLikelihoodMultinomProduct(theSameORF, z) + dlnorm(z,meanlog,sdlog,log=log)
}

unscaledPosterior <- function(theSameORF, z, meanlog, sdlog,log = TRUE){
  exp(logUnscaledPosterior(theSameORF, z,meanlog,sdlog,log = TRUE))
}

#This section is a function exponentiating the above result to attain the product of the multinomial coefficient and the likelihood in its regular (non-logarithmic) form
LikelihoodMultinomProduct <- function(theSameORF, z){
  exp(logLikelihoodMultinomProduct(theSameORF, z))
}
##############################################################################################
#these are functions used to attain the integral of the likelihood function
vec.LnLikelihoodMultinomProduct <- Vectorize(logLikelihoodMultinomProduct,"z")
vec.LikelihoodMultinomProduct <- Vectorize(LikelihoodMultinomProduct,"z")#vectorized list of likelihood values for particular gene with varied z values
vec.unscaledPosterior <- Vectorize(unscaledPosterior,"z")
vec.logUnscaledPosterior <- Vectorize(logUnscaledPosterior,"z")
#This vectorization guarantee that the integration works
#aside from being used in the integrateLikelihood function, vec.getLikelihoodForAnORF can be used to attain the different likelihood values for a specific gene varied with different z value

for ( i in 1:3)
{
  astr<-paste("D:/cygwin/home/Administrator/REU12/scratch/liwe01/FasterVersion/FasterVersion2/test",i,".pdf",sep = "")
  pdf(astr)
  plot(z,vec.LikelihoodMultinomProduct(theSameORFList[[i]],z))
  dev.off()
  
}


plotFunction <- function(z,theSameORF){
  plot(z,vec.LikelihoodMultinomProduct(theSameORF,z))
  
}
#This function integrates the product of the likelihood function and multinomial coefficient
integrateLikelihoodMultinomProduct <- function(theSameORF,lowerb,upperb,rel.tol = .Machine$double.eps^0.5,abs.tol = 0){
  integrand <- function(z){
    vec.LikelihoodMultinomProduct(theSameORF,z)
  }
  integrate(integrand,lowerb,upperb,rel.tol = .Machine$double.eps^0.5,abs.tol = 0)
}

integrateUnscaledPosterior <- function(theSameORF,meanlog,sdlog,log = TRUE,lowerb,upperb,rel.tol = .Machine$double.eps^0.5,abs.tol = 0){
  integrand <- function(z){
    vec.unscaledPosterior(theSameORF, z, meanlog, sdlog,log = TRUE)
  }
  integrate(integrand,lowerb,upperb,rel.tol = .Machine$double.eps^0.5,abs.tol = 0)
}

#Normalizing constant
NormalizingConstantFlatPrior <- function(theSameORF, lowerb, upperb,rel.tol = .Machine$double.eps^0.5,abs.tol = 0){
  1/integrateLikelihoodMultinomProduct(theSameORF,lowerb,upperb,rel.tol = .Machine$double.eps^0.5,abs.tol = 0)$value
}

NormalizingConstant <- function(theSameORF,meanlog,sdlog,log = TRUE,lowerb,upperb,rel.tol = .Machine$double.eps^0.5,abs.tol = 0){
  1/integrateUnscaledPosterior(theSameORF,meanlog,sdlog,log = TRUE,lowerb,upperb,rel.tol = .Machine$double.eps^0.5,abs.tol = 0)$value
}
##########################################################################
#variation of integrateLikelihood function above so that the integral includes the product of the likelihood and z rather than just the likelihood; corresponds to avg z value per each gene
AvgZValuePerGene <- function(theSameORF,lowerb,upperb){
  integrand <- function(z){
    vec.LikelihoodMultinomProduct(theSameORF,z)*z
  }
  integrate(integrand,lowerb,upperb)
}

##########################################################################################################
#This function generates a posterior distribution for a particular gene by dividing the likelihood by the integral of the likelihood(assuming flat prior)

PosteriorExcludingPrior <- function(theSameORF, z,lowerb,upperb){
  integrals <- integrateLikelihoodMultinomProduct(theSameORF,lowerb,upperb)$value
  result <- LikelihoodMultinomProduct(theSameORF, z)
  return(result/integrals)
}

#This is vectorizing the posterior so that it can be integrated
vec.PosteriorExcludingPrior <- Vectorize(PosteriorExcludingPrior,"z")

#This function integrates the likelihood function multiplied by z.  Represents the unscaled mean z value.
unscaledMeanZvalueFlatPrior <- function(theSameORF,lowerb,upperb,rel.tol = .Machine$double.eps^0.5,abs.tol = 0){
  integrand <- function(z){
    vec.LikelihoodMultinomProduct(theSameORF,z)*z
  }
  integrate(integrand,lowerb,upperb,rel.tol = .Machine$double.eps^0.5,abs.tol = 0)
}

unscaledMeanZvalue <- function(theSameORF,lowerb,upperb,meanlog,sdlog,log = TRUE,rel.tol = .Machine$double.eps^0.5,abs.tol = 0){
  integrand <- function(z){
    vec.unscaledPosterior(theSameORF, z, meanlog, sdlog,log = TRUE)*z
  }
  integrate(integrand,lowerb,upperb,rel.tol = .Machine$double.eps^0.5,abs.tol = 0)
}


#This function attains the mean z value for a gene by multiplying the unscaled mean z value by the normalizing constant
MeanZvalueFlatPrior <- function(theSameORF,lowerb,upperb,rel.tol = .Machine$double.eps^0.5,abs.tol = 0){
  NormalizingConstantFlatPrior(theSameORF,lowerb,upperb,rel.tol = .Machine$double.eps^0.5,abs.tol = 0)*unscaledMeanZvalueFlatPrior(theSameORF,lowerb,upperb,rel.tol = .Machine$double.eps^0.5,abs.tol = 0)$value
}

MeanZvalue <- function(theSameORF,lowerb,upperb,meanlog,sdlog,log = TRUE,rel.tol = .Machine$double.eps^0.5,abs.tol = 0){
  NormalizingConstant(theSameORF,meanlog,sdlog,log = TRUE,lowerb,upperb,rel.tol = .Machine$double.eps^0.5,abs.tol = 0)*unscaledMeanZvalue(theSameORF,lowerb,upperb,meanlog,sdlog,log = TRUE,rel.tol = .Machine$double.eps^0.5,abs.tol = 0)$value
}
#This function attains the unscaled mean z-squared value to be used in the variance of z function
UnscaledMeanZSquaredvalueFlatPrior <- function(theSameORF,lowerb,upperb){
  integrand <- function(z){
    vec.LikelihoodMultinomProduct(theSameORF,z)*z^2
  }
  integrate(integrand,lowerb,upperb)
}

#This function attains the mean z squared value for a gene by multiplying the unscaled mean z squared value by the normalizing constant
MeanZsquaredFlatPrior <- function(theSameORF,lowerb,upperb){
  NormalizingConstant(theSameORF,lowerb,upperb)*UnscaledMeanZSquaredvalueFlatPrior(theSameORF,lowerb,upperb)$value
}
#This function calculates the variance of z in a particular gene
Zvariance <- function(theSameORF, lowerb, upperb){
  (MeanZsquaredFlatPrior(theSameORF, lowerb, upperb)) - (MeanZvalueFlatPrior(theSameORF, lowerb, upperb))^2
}

##############################################################################
#This section attains the MLE, or mode, z value for a given gene
### load the workspace
#load("C:/cygwin/home/dillonj3/REU12/scratch/chaij/workspaceforlikelihood.RData")
require(stats)
### genename: gene; start_pt: the point of z to start optimization; lowerb: lowerbound(0)
### upperb: upper bound; trace: print searching step by step if =1, no print if =0
meanlog <- -1.972672
sdlog <- 1.799923
Find.MAP.z.gene <- function(theSameORF,start_pt,lowerb,upperb,trace=0,meanlog,sdlog,log=TRUE){
  op.func <- function(z)
    -vec.logUnscaledPosterior(theSameORF, z,meanlog,sdlog,log)
    ans <- nlminb(start_pt,op.func,lower=lowerb,upper=upperb,control=list(trace=trace))$par
}

Find.MLE.z.gene <- function(theSameORF,start_pt,lowerb,upperb,trace=0){
  op.func <- function(z)
    -vec.LnLikelihoodMultinomProduct(theSameORF, z)
  ans <- nlminb(start_pt,op.func,lower=lowerb,upper=upperb,control=list(trace=trace))$par
}
#generate a phivalue table with increasing phi value
phivalueFile <-"D:/cygwin/home/Administrator/REU12/data/input_files/phi_values/S.cerevisiae.beyer.phi.csv"
phivalueTable <-read.table(phi,sep = ",")
reorderedPhivalueTable <- phivalueTable[with(phivalueTable, order(V2,V1)),]

### sample code to run
#mle.z3 <- Find.mle.z.gene(ORFnamesUnique[3],1e-4,0,100,1) #find mle of z for the 3rd gene

#This section matches the MLE z values to the corresponding phi values for each gene
#zModeRes<-c()
phiValRes<-c()
for (i in 1:4590)
{
 # zModeRes[i] <- MLEres[[i]]$par
  
  phiValRes[i] <- phivalueTable[i,2]
}
simu.MLE.zero <-c()
simu.MLE.nonzero.list <-c()
for (i in 1:length(simu.MLE.nonzero))
{
  simu.MLE.nonzero.list[i] <- simu.MLE.nonzero[[i]]
}

rescaled.simu.MLE <-c()
for ( i in 1:length(simu.MLERes))
{
  if (simu.MLERes[[i]]==0)
  {
    rescaled.simu.MLE[i] <-1e-4
  }
  else
  {
    rescaled.simu.MLE[i]<-simu.MLERes[[i]]
  }
}

simu.MeanZvalueList <-c()
for ( i in 1:length(simu.MeanZvalueRes))
{
  simu.MeanZvalueList[i]<-simu.MeanZvalueRes[[i]]
}