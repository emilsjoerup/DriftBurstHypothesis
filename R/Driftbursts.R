##################
# Written by Emil Sjoerup, a masters of economics student at the university of Aarhus
# Based on code written by Kim Christensen, a professor at the university of Aarhus
# December 2018
# This code calculates the drift-burst test statistic of Christensen et.al (2018) (working paper as of writing this)
# The neccessary inputs are: Tick frequency log prices and time-stamps of these.
# Please e-mail suggestions for improvements to: emilsjoerup@live.dk
##################


drift_bursts = function(time, logprices, testtimes = seq(34200, 57600, 60),
                         PreAverage = 5, AcLag = -1L, Mean_bandwidth = 300L, 
                         Variance_bandwidth = 900L, bParallelize = FALSE, iCores = NA){
  ###Checks###
  if (Mean_bandwidth < 0 | Mean_bandwidth%%1!=0) {
    stop("Mean_bandwidth must be a positive integer")
  }
  if(Variance_bandwidth<0 | Variance_bandwidth%%1!=0){
    stop("Variance_bandwidth must be a positive integer")
  }
  if(AcLag !=-1 && AcLag%%1!=0 || -1>AcLag){
    stop("AcLag must be a positive integer, the standard of -1 designated usage of an automated lag selection algorithm.")
  }
  if(is.na(iCores) & bParallelize){
    print("No iCores argument was provided, sequential evaluation is used.")
  }
  if(anyNA(c(time , logprices , testtimes))){
    stop("NA's in time, logprices or testtimes - might cause crashes and are thus disallowed")
  }
  ###Checks end###
  
  #########  Initialization  ############
  k              = PreAverage #Greatly improves readability at virtually no speed-loss.
  AcLag          = AcLag
  vX             = diff(logprices)
  iT             = length(logprices)
  vPreAveraged   = rep(0, iT-1)
  #########  init end  ############
  
  vPreAveraged[(k*2-1):(iT-1)] = filter(x = logprices, c(rep(1,k),rep(-1,k)))[k:(iT-k)] #Preaveraging
  
  if(bParallelize & !is.na(iCores)){ #Parallel evaluation or not?
   lDriftBursts = DriftBurstLoopCPAR(vPreAveraged, vX, time, testtimes, Mean_bandwidth, Variance_bandwidth, PreAverage, AcLag, iCores )
  }
  else{
   lDriftBursts = DriftBurstLoopC(vPreAveraged, vX, time, testtimes, Mean_bandwidth, Variance_bandwidth, PreAverage, AcLag  )  
  }
  if(is.infinite(lDriftBursts[["DriftBursts"]])[1]){
    cat("unknown error happened in C++ code.")
    return(NULL)
  }
  
  lDriftBursts[["DriftBursts"]][1] = 0 
  #Test cannot be calculated before the first period of testtimes has passed. 
  #Note also, it is unstable until Variance_bandwidth time has passed, but can be calculated.
  
  return(lDriftBursts)
  ## change to return list of mean, variance, and aclags
  
}


