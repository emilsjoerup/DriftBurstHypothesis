driftBursts = function(timestamps = NULL, logPrices, testTimes = seq(34200, 57600, 60),
                        preAverage = 5, ACLag = -1L, meanBandwidth = 300L, 
                        varianceBandwidth = 900L, parallelize = FALSE, nCores = NA, warnings = TRUE){

  #########  Initialization  ############
  k                    = preAverage 
  vX                   = diff(logPrices)
  iT                   = length(logPrices)
  vpreAveraged         = rep(0, iT - 1)
  xts                  = FALSE
  pad = removedFromEnd = 0
  #########  init end  ############
  
  ###Checks###
  if (meanBandwidth<0 | meanBandwidth %% 1 != 0) {
    stop("meanBandwidth must be a positive integer")
  }
  if(varianceBandwidth<0 | varianceBandwidth %% 1 != 0){
    stop("varianceBandwidth must be a positive integer")
  }
  if(ACLag !=-1 && ACLag%%1!=0 | -1>ACLag | ACLag == 1){
    stop("ACLag must be a positive integer greater than or equal to 1, or -1. \n
         The standard of -1 designates usage of an automated lag selection algorithm.")
    #Specifically Newey-West 1994
  }
  if(preAverage <=0 | preAverage%%1!=0 ){
    stop("preAverage must be a positive integer. No preaveraging is done when preAverage = 1.")
  }
  if(inherits(logPrices, "xts")){
    tz        = tzone(logPrices)
    timestamps      = index(logPrices)
    timestamps      = as.numeric(timestamps) - (.indexDate(logPrices)[1] * 86400)
    vIndex    = as.POSIXct((.indexDate(logPrices)[1] * 86400) + testTimes, origin = "1970-01-01")
    logPrices = as.numeric(t(logPrices)) ##need to transpose, otherwise the program will crash.
    vX        = as.numeric(vX)[-1] ### need to remove first entry because diff() on an xts object produces NA in first entry.
    xts       = TRUE
  }
  if((anyNA(timestamps) & !is.null(timestamps)) | anyNA(logPrices) | anyNA(testTimes)){
    stop("NA's in timestamps, logPrices or testTimes - might cause crashes and are thus disallowed")
  }
  if((length(timestamps) != length(logPrices) & !is.null(timestamps))){
    stop("timestamps and logPrices input not of same length, to prevent crashing this is not allowed.")
  }
  if((is.na(nCores) | nCores %% 1 != 0) & parallelize){
    warning("No iCores argument was provided, or the provided nCores argument is not an integer.\n
          Sequential evaluation is used.")
    bParallelize = FALSE
  }
  if(15 <= max(diff(timestamps))/60){
    stop("There is a period of greater than 15 minutes with no trades.\n
         In my testing this may cause crashes and is thus disallowed")
  }
  if(min(timestamps) > min(testTimes[-1])){
    testTimes = testTimes[-2]
    pad = 1
    while(min(timestamps) > min(testTimes[-1])){
      testTimes = testTimes[-2] 
      pad = pad + 1
    }
    if(warnings){
    warning(paste("\nThe first testing timestamps is  before any observations. May cause crashes.",
                  "\nItereatively removing first testing timestamps until this is no longer the case.",
                  "\nremoved the first", pad, "entries from testTimes and replacing with a 0\n"))
    }
  }
  if(max(testTimes)>max(timestamps) + 900){
    testTimes = testTimes[-length(testTimes)]
    removedFromEnd = 1
    while(max(testTimes)>max(timestamps) + 900){
      testTimes = testTimes[-length(testTimes)]  
      removedFromEnd = removedFromEnd + 1
    }
    if(warnings){
      warning(paste("\nThe last testing timestamps is more than 15 minutes after the last trade, this may cause crashes.",
                    "\nIteratively removing the last testing timestamps until this is no longer the case.",
                    "\nremoved the last", removedFromEnd, "entries from testTimes\n"))
      
    }
    
  }
  ###Checks end###
  

  vpreAveraged[(k*2 - 1):(iT - 1)] = filter(x = logPrices, c(rep(1,k),rep( -1,k)))[k:(iT - k)] #Preaveraging
  
  if(parallelize & !is.na(nCores)){ #Parallel evaluation or not?
   lDriftBursts = DriftBurstLoopCPAR(c(0,vpreAveraged), c(0,vX), timestamps, testTimes, meanBandwidth, 
                                     varianceBandwidth, preAverage, ACLag, nCores )
  }
  else{
   lDriftBursts = DriftBurstLoopC(c(0,vpreAveraged), c(0,vX), timestamps, testTimes, meanBandwidth, 
                                  varianceBandwidth, preAverage, ACLag)  
  }
  
  lDriftBursts[["DriftBursts"]][1] = 0 
  lDriftBursts[["Sigma"]][1]       = 0
  lDriftBursts[["Mu"]][1]          = 0
  
  if(pad != 0 | removedFromEnd != 0){
    lDriftBursts[["DriftBursts"]] = c(rep(0,pad), lDriftBursts[["DriftBursts"]], rep(0,removedFromEnd))
    lDriftBursts[["Sigma"]]       = c(rep(0,pad), lDriftBursts[["Sigma"]], rep(0,removedFromEnd))
    lDriftBursts[["Mu"]]          = c(rep(0,pad), lDriftBursts[["Mu"]], rep(0,removedFromEnd))
  }
  
  if(xts){
    lDriftBursts[["DriftBursts"]] = xts(lDriftBursts[["DriftBursts"]], order.by = vIndex, tz = tz)
    lDriftBursts[["Sigma"]]       = xts(lDriftBursts[["Sigma"]], order.by = vIndex, tz = tz)
    lDriftBursts[["Mu"]]          = xts(lDriftBursts[["Mu"]], order.by = vIndex, tz = tz)
  }
  
  lInfo = list("varianceBandwidth" = varianceBandwidth, "meanBandwidth" = meanBandwidth,"preAverage" = preAverage,
               "nObs" = iT, "testTimes" = testTimes)
  lDriftBursts[["Info"]] = lInfo
  #replace NANs with 0's
  NANS = is.nan(lDriftBursts[["Sigma"]])
  lDriftBursts[["DriftBursts"]][NANS] = 0
  lDriftBursts[["Sigma"]][NANS]       = 0
  
  class(lDriftBursts) = c("DBH", "list")
  return(lDriftBursts)
}

drift_bursts = function(time = NULL, logprices, testtimes = seq(34200, 57600, 60),
                        PreAverage = 5, AcLag = -1L, Mean_bandwidth = 300L, 
                        Variance_bandwidth = 900L, bParallelize = FALSE, iCores = NA, warnings = TRUE){
  .Deprecated("driftBursts", msg = "This package has moved to camelCase, drift_bursts have been renamed to driftBursts and arguments have been renamed")
  return(driftBursts(time, logprices, testtimes, PreAverage, AcLag, Mean_bandwidth, Variance_bandwidth,
                     bParallelize, iCores, warnings))
}