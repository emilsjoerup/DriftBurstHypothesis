##################
# Written by Emil Sjoerup, a masters of economics student at the university of Aarhus
# Based on code written by Kim Christensen, a professor at the university of Aarhus
# December 2018
# This code calculates the drift-burst test statistic of Christensen et.al (2018) 
# (working paper as of writing this)
# The neccessary inputs are: Tick frequency log prices and time-stamps of these.
# Please e-mail suggestions for improvements to: emilsjoerup@live.dk
##################

drift_bursts = function(time = NULL, logprices, testtimes = seq(34200, 57600, 60),
                        PreAverage = 5, AcLag = -1L, Mean_bandwidth = 300L, 
                        Variance_bandwidth = 900L, bParallelize = FALSE, iCores = NA, warnings = TRUE){

  #########  Initialization  ############
  k              = PreAverage 
  AcLag          = AcLag
  vX             = diff(logprices)
  iT             = length(logprices)
  vPreAveraged   = rep(0, iT-1)
  xts            = FALSE
  pad            = 0
  #########  init end  ############
  
  ###Checks###
  if (Mean_bandwidth<0 | Mean_bandwidth%%1!=0) {
    stop("Mean_bandwidth must be a positive integer")
  }
  if(Variance_bandwidth<0 | Variance_bandwidth%%1!=0){
    stop("Variance_bandwidth must be a positive integer")
  }
  if(AcLag !=-1 && AcLag%%1!=0 | -1>AcLag){
    stop("AcLag must be a positive integer or -1, the standard of -1 designates usage of an automated lag selection algorithm.")
    #Specifically Newey-West 1994
  }
  if(inherits(logprices, "xts")){
    time = index(logprices)
    time = as.numeric(time) - (.indexDate(logprices)[1] * 86400)
    vIndex = as.POSIXct((.indexDate(logprices)[1] * 86400) + testtimes, origin = "1970-01-01")
    logprices = as.numeric(t(logprices)) ##need to transpose, otherwise the program will crash.
    vX = as.numeric(vX)[-1] ### need to remove first entry because diff() on an xts object produces NA in first entry.
    xts = TRUE
  }
  if((anyNA(time) & !is.null(time)) | anyNA(logprices) | anyNA(testtimes)){
    stop("NA's in time, logprices or testtimes - might cause crashes and are thus disallowed")
  }
  if((length(time) != length(logprices) & !is.null(time))){
    stop("Time and logprices input not of same length, to prevent crashing this is not allowed.")
  }
  if((is.na(iCores) | iCores%%1!=0) & bParallelize ){
    print("No iCores argument was provided, or the provided iCores argument is not an integer. Sequential evaluation is used.")
    bParallelize = FALSE
  }
  if(min(time) > min(testtimes[-1])){
    testtimes = testtimes[-2]
    pad = 1
    while(min(time) > min(testtimes[-1])){
      testtimes = testtimes[-2] 
      pad = pad+1
    }
    if(warnings){
    cat("\nThe first testing time is  before any observations. May cause crashes, thus disallowed.")
    cat("\nItereatively removing first testing time until this is no longer the case.")
    cat("removed", pad, "entries from testtimes")
    }
  }
  
  
  ###Checks end###
  vPreAveraged[(k*2-1):(iT-1)] = filter(x = logprices, c(rep(1,k),rep(-1,k)))[k:(iT-k)] #Preaveraging

  if(bParallelize & !is.na(iCores)){ #Parallel evaluation or not?
   lDriftBursts = DriftBurstLoopCPAR(vPreAveraged, vX, time, testtimes, Mean_bandwidth, 
                                     Variance_bandwidth, PreAverage, AcLag, iCores )
  }
  else{
   lDriftBursts = DriftBurstLoopC(vPreAveraged, vX, time, testtimes, Mean_bandwidth, 
                                  Variance_bandwidth, PreAverage, AcLag)  
  }
  
  lDriftBursts[["DriftBursts"]][1] = 0 
  lDriftBursts[["Sigma"]][1]       = 0
  lDriftBursts[["Mu"]][1]          = 0
  
  if(pad!=0){
    lDriftBursts[["DriftBursts"]] = c(rep(0,pad), lDriftBursts[["DriftBursts"]])
    lDriftBursts[["Sigma"]]       = c(rep(0,pad), lDriftBursts[["Sigma"]])
    lDriftBursts[["Mu"]]          = c(rep(0,pad), lDriftBursts[["Mu"]])
  }
  
  if(xts){
    lDriftBursts[["DriftBursts"]] = xts(lDriftBursts[["DriftBursts"]], order.by = vIndex )
    lDriftBursts[["Sigma"]]       = xts(lDriftBursts[["Sigma"]], order.by = vIndex )
    lDriftBursts[["Mu"]]          = xts(lDriftBursts[["Mu"]], order.by = vIndex )
  }
  #Test cannot be calculated before the first period of testtimes has passed. 

  return(lDriftBursts)
  
  
}





