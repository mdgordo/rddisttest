#' @title Tests for Manipulation of Discrete Running Variable
#' @description Provides an alternative to conventional RD manipulation tests, such as the McCrary test, to be applied in the case where the running variable takes on discrete values.
#' @param RV Vector containing the values for each observations level of the running variable.
#' @param C Real number specifying the level of the RV at which the discontinuity exists.
#' @param K Real number corresponding to the user specified bound on the RV’s probability mass function curvature.
#'
#' @return P-value
#' @export
#' @importFrom stats pbinom qbinom
#' @examples
#' \dontrun{
#' #rddisttest(dataframe$runningvar, 0, .01)
#' }

rddisttest<-function(RV, C, W, K=0){
  freqdf<- aggregate(x = list("Freq" = W), by = list("RV" = RV), FUN = sum)
  freqdf$rank<-1:nrow(freqdf)
  freqdf$Freq <- round(freqdf$Freq, 0)
  c.rank<-as.numeric(min(freqdf[which(as.numeric(freqdf$RV)>=C),3]))

  #Determine frequencies at, below, and above the threshold
  nc<-as.numeric(freqdf[which(freqdf$rank==c.rank),2])
  ncb<-as.numeric(freqdf[which(freqdf$rank==(c.rank-1)),2])
  nca<-as.numeric(freqdf[which(freqdf$rank==(c.rank+1)),2])
  n<-(nc + ncb + nca)

  #Now, we create a function that tests the null hypothesis at various significance levels.
  #We do this by obtaining a ?X4 matrix containing all potential critical values that
  #meet the constraint of sufficient coverage. After obtaining all candidate pairs, we
  #choose as our critical values the one with maximum coverage amongst those that minimize |Cu-Cl|

  TestH0 <- function(RV, C, K, A) {
    keep<-c(0,n,n,0)
    if (K==0){
      Cu<-n+1
      Cl<-0
      areaCu.new<-1

      #Run the while loop only for upper CVs where area below has sufficient coverage
      while (round(areaCu.new, digits = 10) >= (1-A)){
        areaCu.old<-pbinom(Cu-1,n,1/3)
        Cl<-(qbinom((A+areaCu.old-1),n,1/3)-1)
        areaCl<-pbinom(Cl,n,1/3)
        coverage<-(areaCu.old-areaCl)
        #If coverage is sufficient we keep the candidate pair and check for sufficient coverage at a higher Cl
        if (coverage>=(1-A)){
          addkeep<-c(Cl,Cu,(Cu-Cl),coverage)
          keep<-rbind(keep,addkeep)
          checkcoverage<-areaCu.old - pbinom(Cl+1,n,1/3)
          if (checkcoverage >= (1-A)){
            addkeep<-c(Cl+1,Cu,Cu-Cl-1,checkcoverage)
            keep<-rbind(keep,addkeep)
          }
        }
        #In the case where coverage is not sufficient we check to see if it is sufficient when we drop Cl by 1
        else{
          checkcoverage<-areaCu.old - pbinom(Cl-1,n,1/3)
          if (checkcoverage >= (1-A)) {
            addkeep<-c(Cl-1,Cu,Cu-Cl+1,checkcoverage)
            keep<-rbind(keep,addkeep)
          }
        }
        Cu<-Cu-1
        areaCu.new<-pbinom(Cu-1,n,1/3)
      }
    }

    #The process is similar, yet slightly more complicated when k!=0. The main difference being
    #that now we only consider candidate pairs that satisfy the minimum coverage constraint at
    #both ends of the success probability interval
    else {
      plow<-((1-K)/(3-K))
      phigh<-((1+K)/(3+K))
      areaCu.newL<-1
      areaCu.newH<-1
      Cu<-n+1
      while (round(areaCu.newL, digits = 10)>=(1-A) & round(areaCu.newH, digits = 10)>=(1-A)){
        areaCu.oldL<-pbinom(Cu-1,n,plow)
        areaCu.oldH<-pbinom(Cu-1,n,phigh)
        Cl<-(min(c(qbinom((A+areaCu.oldL-1),n,plow),qbinom((A+areaCu.oldH-1),n,phigh)))-1)
        coverage.L<-(areaCu.oldL-pbinom(Cl,n,plow))
        coverage.H<-(areaCu.oldH-pbinom(Cl,n,phigh))
        coverage<-min(coverage.L,coverage.H)
        if (coverage >= (1-A)){
          addkeep<-c(Cl,Cu,Cu-Cl,coverage)
          keep<-rbind(keep,addkeep)
          check.covL<-areaCu.oldL - pbinom(Cl+1,n,plow)
          check.covH<-areaCu.oldH - pbinom(Cl+1,n,phigh)
          checkcoverage<-min(check.covL,check.covH)
          if (checkcoverage >= (1-A)){
            addkeep<-c(Cl+1,Cu,Cu-Cl-1,checkcoverage)
            keep<-rbind(keep,addkeep)
          }
        }
        else {
          check.covL<-areaCu.oldL - pbinom(Cl-1,n,plow)
          check.covH<-areaCu.oldH - pbinom(Cl-1,n,phigh)
          checkcoverage<-min(check.covL,check.covH)
          if (checkcoverage>=(1-A)){
            addkeep<-c(Cl-1,Cu,Cu-Cl+1,checkcoverage)
            keep<-rbind(keep,addkeep)
          }
        }

        Cu<-Cu-1
        areaCu.newL<-pbinom(Cu-1,n,plow)
        areaCu.newH<-pbinom(Cu-1,n,phigh)


      }
    }
    leastdif<-keep[which.min(keep[,3]),3]
    Minimizers<-keep[which(keep[,3]==leastdif),]
    #If the minimum is unique, we have our critical values
    if (NCOL(Minimizers)==1){
      Winner<-Minimizers
    }
    #If the minimum is not unique, then we choose our critical values to be those that maximize coverage.
    else{
      Winner<-Minimizers[which.max(Minimizers[,4]),]
    }
    Cl.star<-Winner[1]
    Cu.star<-Winner[2]
    if (nc<Cu.star & nc>Cl.star){
      reject<-0
    }
    else {
      reject<-1
    }
    invisible(reject)
  }

  #Use a bisecting procedure for ten iterations to approximate the p-value
  L<-1:10
  A<-.5
  delta<-.5
  for (i in L){
    if (TestH0(RV,C,K,A)==1){
      delta<-delta/2
      A<-A-(delta)

    }
    else {
      delta<-delta/2
      A<-A+(delta)
    }
  }
  p.round<-round(A, digits=3)


  cat("RD Density Test Results, K =",K,"\n","p-value =",p.round)
  invisible(A)
}


