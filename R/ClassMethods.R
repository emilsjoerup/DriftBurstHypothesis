plot.DBH = function(x, ...){
  
  #### Get extra passed options and data
  options = list(...)
  #### List of standard options
  opt = list(which = "driftbursts", price = NULL, timestamps = NULL, startTime = 34200, endTime = 57600,
             leg.x = "topleft", leg.y = NULL,  tz = "GMT", annualize = FALSE, nDays = 252, legend.txt = "")
  #### Override standard options where user passed new options
  opt[names(options)] = options
  # for (i in 1:length(opt)) {
  #   assign(names(opt)[i], opt[[i]])
  # }

  #### Extract options (better way to do it?)  #I will test above method when I get time.
  which      = tolower(opt$which)
  startTime  = opt$startTime
  endTime    = opt$endTime
  main       = opt$main
  tz         = opt$tz
  leg.x      = opt$leg.x
  leg.y      = opt$leg.y
  annualize  = opt$annualize
  nDays      = opt$nDays
  price      = opt$price
  timestamps = opt$timestamps
  tstat      = getDB(x)
  sigma      = getSigma(x, annualize, nDays) 
  mu         = getMu(x, annualize, nDays)
  MB         = x$Info$meanBandwidth
  VB         = x$Info$varianceBandwidth
  startpar   = par(no.readonly = TRUE)
  testTimes  = x$Info$testTimes
  horizLines = seq(round(min(tstat)), round(max(tstat)), 1)
  ###Setup done
  if(!all(which %in% c("driftbursts", "mu", "sigma", "db"))){
    stop("The which argument must be a character vector containing either:\n
         Sigma, Mu, both of these or DriftBursts. 
         CasE doesn't matter.")
  }
  if(inherits(tstat, "xts")){
    testTimes = index(tstat)
    testTimes = as.numeric(testTimes) - (.indexDate(tstat)[1] * 86400)
    tstat     = as.numeric(tstat)
    sigma     = as.numeric(sigma)
    mu        = as.numeric(mu)
    if(!is.null(price)){
      timestamps  = index(price)
      timestamps  = as.numeric(timestamps) - (.indexDate(price)[1] * 86400)
      price       = as.numeric(price)
    }
  }
  if(testTimes[1] == startTime){
    testTimes = testTimes[-1]
    sigma = sigma[-1]
    mu = mu[-1]
    tstat=tstat[-1]
  }
  xtext = as.POSIXct(testTimes, origin = "1970-01-01", tz = tz)
  xlim = c(startTime, endTime)
  xlab = "Time"
  
  if(all(which %in% c("driftbursts", "db"))){ #use all() because this function should accept which arguments with length longer than 1
    par(mar = c(4,3.5,2,1.25), mgp = c(2,1,0))
    if(!is.null(price)) par(mar = c(4,3.5,4,4), mgp = c(2,1,0)) #makes room for values on the right y-axis
    main = "Drift Bursts test statistic"
    ylab = "test-statistic"
    plot(tstat, x = xtext, type = "l", xaxt = 'n', ylab = ylab, main = main, xlab = xlab, xlim = xlim)
    axis.POSIXct(side  = 1, at = seq.POSIXt(xtext[1], xtext[length(xtext)], length.out = 7))
    abline(h = horizLines, col = "grey" , lty = 3, cex = 0.1)
    legend.txt = "t-stat"
    if(!is.null(price)){
      if(is.null(timestamps)){
        stop("The timestamps of the price must be passed in the timestamps argument")
      }
      par(new = TRUE)
      plot(price, x = timestamps , type = "l", axes = FALSE, col = "red", xlab = "", ylab = "", lty = 2, xlim = xlim)  
      axis(4)
      mtext(side = 4, text = "price", line = 2.5)
      legend.txt = c(legend.txt, "price")
      legend(x = leg.x, leg.y, legend = legend.txt, lty = c(1,2), col = c(1,2), bg = rgb(0,0,0,0), box.lwd = 0, 
             box.col = rgb(0,0,0,0))
    }
  }
  if(all(which == "sigma")){ #use all() because this function should accept which arguments with length longer than 1
    main = "volatility"
    ylab = "local volatility"
    par(mar = c(4,3.5,2,1.25), mgp = c(2,1,0))
    plot(sigma, x = xtext, type = "l",  xaxt = 'n', ylab = ylab, main = main, xlab = xlab)  
    axis.POSIXct(side  = 1, at = seq.POSIXt(xtext[1], xtext[length(xtext)], length.out = 7))
  }
  if(all(which == "mu")){ #use all() because this function should accept which arguments with length longer than 1
    main = "drift"
    ylab = "drift"
    if(annualize){ 
      ylab = paste("annualized", ylab)
    }
    par(mar = c(4,3.5,2,1.25), mgp = c(2,1,0))
    plot(mu, x = xtext, type = "l",  xaxt = 'n', ylab = ylab, main = main, xlab = xlab)
    axis.POSIXct(side  = 1, at = seq.POSIXt(xtext[1], xtext[length(xtext)], length.out = 7))
    abline(h = 0, col = "grey" , lty = 3)
  }
  if("mu" %in% which & "sigma" %in% which){
    par(mfrow = c(2,1), omi = c(0,0,0,0), mgp = c(2,1,0), mai = c(0.75,0.75,0.3,0.25))
    main = "drift"
    ylab = "drift"
    if(annualize){ 
      ylab = paste("annualized", ylab)
    }
    plot(mu, x = xtext, type = "l", xlab = "",  xaxt = 'n', ylab = ylab, main = main)
    axis.POSIXct(side  = 1, at = seq.POSIXt(xtext[1], xtext[length(xtext)], length.out = 7))
    abline(h = 0, col = "grey" , lty = 3)
    main = "volatility"
    ylab = "volatility"
    if(annualize){ 
      ylab = paste("annualized", ylab)
    }
    plot(sigma, x = xtext, type = "l", xlab = "", xaxt = 'n', ylab = ylab, main = main) 
    axis.POSIXct(side  = 1, at = seq.POSIXt(xtext[1], xtext[length(xtext)], length.out = 7))
  }
  par(startpar)
}

print.DBH = function(x, ...){
  options = list(...)
  #### List of standard options
  opt = list(annualize = FALSE, nDays = 252)
  #### Override standard options where user passed new options
  opt[names(options)] = options
  
  annualize = opt$annualize
  nDays     = opt$nDays
  
  cat("\nMean mu:                   ", mean(getMu(x, annualize, nDays)[-1]))
  cat("\nMean sigma:                ", mean(getSigma(x, annualize, nDays)[-1]))
  cat("\nMean test statistic:       ", mean(getDB(x)[-1]))
  cat("\nVariance of test statistic:", var(getDB(x)[-1]))
}

getDB = function(object){
  UseMethod("getDB", object)
}

getDB.DBH = function(object){
  DB = object$driftBursts
  return(DB)
}


getSigma = function(object, annualize = FALSE, nDays = 252){
  UseMethod("getSigma", object)
}

getSigma.DBH = function(object, annualize = FALSE, nDays = 252){
  sigma = sqrt((object$sigma * 2 * object$Info$nObs)  / (object$Info$nObs / 23400))/(object$Info$preAverage^2)
  if(annualize){sigma = sigma * sqrt(nDays)}
  return(sigma)
}


getMu = function(object, annualize = FALSE, nDays = 252){
  UseMethod("getMu", object)
}

getMu.DBH = function(object, annualize = FALSE, nDays = 252){
  mu = (object$mu * object$Info$meanBandwidth / (object$Info$nObs / 23400)) / (object$Info$preAverage^2 * 2)
  if(annualize){mu =  mu * nDays}
  return(mu)
}

