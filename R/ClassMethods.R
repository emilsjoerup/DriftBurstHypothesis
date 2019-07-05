plot.DBH = function(x, ...){
  
  #### Get extra passed options and data
  options = list(...)
  #### List of standard options
  opt = list(which = "DriftBursts", price = NULL, time = NULL, startTime = 34200, endTime = 57600,
  xlab = "time", ylab = "value",  leg.x = "topleft", leg.y = NULL, tz = "GMT", annualize = FALSE, nDays = 252)
  #### Override standard options where user passed new options
  opt[names(options)] = options
  #### Extract options (better way to do it?)
  which = tolower(opt$which)
  price = opt$price
  time  = opt$time
  startTime = opt$startTime  
  endTime = opt$endTime  
  xlab = opt$xlab
  ylab = opt$ylab
  main = opt$main
  tz    = opt$tz
  leg.x = opt$leg.x
  leg.y = opt$leg.y
  annualize = opt$annualize
  nDays     = opt$nDays
  tstat = x$DriftBursts
  sigma = getSigma(x, annualize, nDays) 
  mu    = getMu(x, annualize, nDays)
  MB    = x$Info$Mean_Bandwidth
  VB    = x$Info$Variance_Bandwidth
  startpar = par(no.readonly = TRUE)
  legend.txt = ""
  testtimes  = seq(startTime, endTime, length.out = length(tstat))
  horizLines = seq(round(min(tstat)), round(max(tstat)), 2)
  ###Setup done
  if(!all(which %in% c("driftbursts", "mu", "sigma"))){
    stop("The which argument must be a character vector containing either:\n
         Sigma, Mu, both of these or DriftBursts. 
         Case doesn't matter.")
  }
  if(inherits(tstat, "xts")){
    testtimes = index(tstat)
    testtimes = as.numeric(testtimes) - (.indexDate(tstat)[1] * 86400)
    tstat = as.numeric(tstat)
    sigma = as.numeric(sigma)
    mu    = as.numeric(mu)
    if(!is.null(price)){
      time = index(price)
      time = as.numeric(time) - (.indexDate(price)[1] * 86400)
      price = as.numeric(price)
    }
  }
  xtext = as.POSIXct(testtimes, origin = "1970-01-01", tz = tz)
  if(all(which %in% c("driftbursts", "db"))){ #use all() because this function should accept which arguments with length longer than 1
    par(mar = c(4,3.5,2,1.25), mgp = c(2,1,0))
    if(!is.null(price)) par(mar = c(4,3.5,4,4), mgp = c(2,1,0)) #makes room for values on the right y-axis
    main = "Drift Bursts test statistic"
    plot(tstat, x = xtext, type = "l", xaxt = 'n', ylab = ylab, main = main, xlab = xlab)
    axis.POSIXct(side  = 1, at = seq.POSIXt(xtext[1], xtext[length(xtext)], length.out = 7))
    abline(h = horizLines, col = "grey" , lty = 3)
    legend.txt = "t-stat"
    if(!is.null(price)){
      if(is.null(time)){
        stop("The timestamps of the price must be passed in the time argument")
      }
      par(new = TRUE)
      plot(price, x = time/86400 , type = "l", axes = FALSE, col = "red", xlab = "", ylab = "", lty = 2)  
      axis(4)
      mtext(side = 4, text = "price", line = 2.5)
      legend.txt = c(legend.txt, "price")
      legend(x = leg.x, leg.y, legend = legend.txt, lty = c(1,2), col = c(1,2), bg = rgb(0,0,0,0), box.lwd = 0, 
           box.col = rgb(0,0,0,0))
    }
  }
  if(all(which == "sigma")){ #use all() because this function should accept which arguments with length longer than 1
    main = "Sigma"
    par(mar = c(4,3.5,2,1.25), mgp = c(2,1,0))
    plot(sigma, x = xtext, type = "l",  xaxt = 'n', ylab = ylab, main = main, xlab = xlab)  
    axis.POSIXct(side  = 1, at = seq.POSIXt(xtext[1], xtext[length(xtext)], length.out = 7))
  }
  if(all(which == "mu")){ #use all() because this function should accept which arguments with length longer than 1
    main = "Mu"
    par(mar = c(4,3.5,2,1.25), mgp = c(2,1,0))
    plot(mu, x = xtext, type = "l",  xaxt = 'n', ylab = ylab, main = main, xlab = xlab)
    axis.POSIXct(side  = 1, at = seq.POSIXt(xtext[1], xtext[length(xtext)], length.out = 7))
    abline(h = 0, col = "grey" , lty = 3)
  }
  if("mu" %in% which & "sigma" %in% which){
    par(mfrow = c(2,1), omi = c(0,0,0,0), mgp = c(2,1,0), mai = c(0.75,0.75,0.3,0.25))
    main = "Mu"
    plot(mu, x = xtext, type = "l", xlab = "",  xaxt = 'n', ylab = ylab, main = main)
    axis.POSIXct(side  = 1, at = seq.POSIXt(xtext[1], xtext[length(xtext)], length.out = 7))
    abline(h = 0, col = "grey" , lty = 3)
    main = "Sigma"
    plot(sigma, x = xtext, type = "l", xlab = "", xaxt = 'n', ylab = ylab, main = main) 
    axis.POSIXct(side  = 1, at = seq.POSIXt(xtext[1], xtext[length(xtext)], length.out = 7))
  }
  par(startpar)
  invisible(NULL)
}



getDB = function(object){
  UseMethod("getDB", object)
}

getDB.DBH = function(object){
  DB = object$DriftBursts
  return(DB)
}


getSigma = function(object, annualize = FALSE, nDays = 252){
  UseMethod("getSigma", object)
}

getSigma.DBH = function(object, annualize = FALSE, nDays = 252){
  sigma = object$Sigma * 2 * object$Info$nObs
  if(annualize){sigma = sigma * sqrt(nDays)}
  return(sigma)
}


getMu = function(object, annualize = FALSE, nDays = 252){
  UseMethod("getMu", object)
}

getMu.DBH = function(object, annualize = FALSE, nDays = 252){
  mu = object$Mu * object$Info$Mean_bandwidth
  if(annualize){mu =  mu * sqrt(nDays)}
  return(mu)
}

