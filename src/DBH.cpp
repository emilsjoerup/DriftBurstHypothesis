#include <RcppArmadillo.h>
#if defined(_OPENMP)
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif
using namespace arma;

const int iMaxLag = 20;

// [[Rcpp::export]]
arma::vec HACWeightC(int iLag){ 
  // Input will always be non-zero and positive, thus, the 0 case is ignored.
  vec vW = linspace(1, iLag, iLag) / iLag;
  int iIdx = sum(vW<=0.5);
  vW(span(0,iIdx-1)) = 1-6*square(vW(span(0,iIdx-1))) + 6 * pow(vW(span(0,iIdx-1)), 3);
  vW(span(iIdx, iLag-1)) = 2 * pow(1 - vW(span(iIdx, iLag-1)), 3);
  return(vW);
}


//[[Rcpp::export]]
double AsymptoticVarianceC(arma::vec vIn, int iLag){
    try{
      double d0 = sum(vIn % vIn);
      vec vAc = zeros(iLag);
      if(iLag == 0){
        vec vAc = zeros(1);
      }
      int iT = vIn.size();
      if(iT<=iLag){
        return(datum::nan);
      }
      for(int i=0; i<iLag; i++){
        vAc[i] = sum(vIn(span((i+1), iT-1)) % vIn(span(0, (iT-i-2))));
      }
      vec vW = HACWeightC(iLag);
      double dOut = d0 + 2.0 * sum(vAc % vW);
      return(dOut);
    }
    catch( std::exception &ex ) {
      Rf_warning("unkown C++ error occured in asymptotic variance\n");
      return(datum::nan);
    } catch(...) {
      Rf_warning("unkown C++ error occured in asymptotic variance\n");
      return(datum::nan);
    }
}


// [[Rcpp::export]]
double AutomaticLagSelectionC(arma::vec vX, double dMu){
  double dC = 2.6614;
  int iN = round(4.0 * pow((dMu/100), 0.16));
  int iT = vX.size();
  if(iT<=iN+1){
    return(iMaxLag);
  }
  vec vE = zeros(iN+1);
  vec vAc = zeros(iN);
  
  for(int i=0; i<(iN+1); i++){
    vE[i] = mean(vX(span((i+1), iT-1)) % vX(span(0, (iT-i-2))));
  }
  
  double d0 = - sum(linspace(1, (iN+1), (iN+1)) % vE);
  
  for(int i=1; i<(iN+1); i++){
    vAc[(i-1)] = sum(linspace(1, (iN-i+1), (iN-i+1)) % vE(span(i, iN)));
  }
  
  double s0 = d0 + 2.0 * -sum(vAc);
  
  
  double sQ = 2.0 * sum(pow(linspace(1, iN, iN), 2.0) % vAc);
  double dGamma = dC * pow(pow((sQ / s0), 2.0 ), 0.2);
  
  int iOut = round(dGamma * pow(dMu, 0.2));
  
  return(iOut);
  
  
}

//[[Rcpp::export]]

Rcpp::List DriftBurstLoopC(arma::vec vPreAveraged, arma::vec diffedlogprices, arma::vec vTime, arma::vec vTesttime, 
                           double iMeanBandwidth, double iVarBandwidth, int iPreAverage, int iAcLag){
  int iT = vTesttime.size();
  int iN = vPreAveraged.size();
  int iQ;
  int iIdx;
  double iAutoAcLag;
  arma::vec vMu = zeros<vec>(iT);
  arma::vec vSigma = zeros<vec>(iT);
  arma::vec vX;
  arma::vec vWm;
  arma::vec vWvar;
  arma::vec vFoo(2);
  vFoo[1] = iMaxLag;
  for(int i=1; i<iT; i++){
    
    vX = vTime(span(0, (iN-1))) - vTesttime[i];
    iIdx = sum((vX<=0))-1;
    
    vWm = exp(-abs(vX(span(0,iIdx))/iMeanBandwidth));  // exponential kernel
    
    vMu[i] = sum(vWm(span(0,iIdx)) % vPreAveraged(span(0,iIdx))) / iMeanBandwidth; 
    
    vWvar = exp(-abs(vX(span(0,iIdx))/iVarBandwidth));  // exponential kernel
    
    if(iAcLag == -1){
      iQ = AutomaticLagSelectionC(diffedlogprices, sum(vWm));
      vFoo[0] = iQ;
      iAutoAcLag = min(vFoo) + 2.0 * (iPreAverage-1);
      vSigma[i] = AsymptoticVarianceC((vWvar(span(0,iIdx)) % vPreAveraged(span(0,iIdx))), iAutoAcLag) / iVarBandwidth; 
      
      
    }else{
      vSigma[i] = AsymptoticVarianceC((vWvar(span(0,iIdx)) % vPreAveraged(span(0,iIdx))), iAcLag) / iVarBandwidth;
  
    }
    
    
  }
  
  arma::vec vDb = sqrt(iMeanBandwidth) * vMu / sqrt(vSigma) ; 
  return(Rcpp::List::create(Rcpp::Named("DriftBursts") = vDb,
                            Rcpp::Named("Sigma") = vSigma,
                            Rcpp::Named("Mu") = vMu));
}




//[[Rcpp::export]]


Rcpp::List DriftBurstLoopCPAR(arma::vec vPreAveraged, arma::vec diffedlogprices, arma::vec vTime, 
                              arma::vec vTesttime, 
                              double iMeanBandwidth, double iVarBandwidth, int iPreAverage, int iAcLag , int iCores){
  
  
#if defined(_OPENMP)
  omp_set_num_threads(iCores);
  int iT = vTesttime.size();
  int iN = vPreAveraged.size();
  int iIdx;
  int i;
  arma::vec vMu= zeros<vec>(iT);
  arma::vec vSigma= zeros<vec>(iT);
  arma::vec vX;
  arma::vec vWm;
  arma::vec vWvar;
  //Parallelization setup
  #pragma omp parallel for default(none)\
  shared(vPreAveraged, diffedlogprices, vTime, vTesttime, iMeanBandwidth,\
         iVarBandwidth, iPreAverage, iAcLag, vMu, vSigma, iT, iN)\
    private(vX, vWm, vWvar, i, iIdx)
    
    for(i=1; i<iT; i++){
      
      vX = vTime(span(0, (iN-1))) - vTesttime[i];
      
      iIdx = sum((vX<=0))-1;
      
      vWm = exp(-abs(vX(span(0,iIdx))/iMeanBandwidth));  // exponential kernel
      
      vMu[i] = sum(vWm(span(0,iIdx)) % vPreAveraged(span(0,iIdx))) / iMeanBandwidth; 
      
      vWvar = exp(-abs(vX(span(0,iIdx))/iVarBandwidth));  // exponential kernel
      
      if(iAcLag == -1){
        
        double iQ = AutomaticLagSelectionC(diffedlogprices, sum(vWm));
        arma::vec vFoo = {iQ, iMaxLag};
        int iAutoAcLag = min(vFoo)+ 2.0 * (iPreAverage-1);
        vSigma[i]= AsymptoticVarianceC((vWvar(span(0,iIdx)) % vPreAveraged(span(0,iIdx))),iAutoAcLag) / iVarBandwidth; 
        
        
      }else{
        vSigma[i] = AsymptoticVarianceC((vWvar(span(0,iIdx)) % vPreAveraged(span(0,iIdx))),iAcLag) / iVarBandwidth;
      }
      
      
      
    }
    
  arma::vec vDb = sqrt(iMeanBandwidth) * vMu / sqrt(vSigma) ; 
  Rcpp::List lOut =  Rcpp::List::create(Rcpp::Named("DriftBursts") = vDb,
                                        Rcpp::Named("Sigma") = vSigma,
                                        Rcpp::Named("Mu") = vMu);
#else
  Rf_warning("OpenMP is not available. Sequential processing is used.");
  lOut = DriftBurstLoopC(vPreAveraged, diffedlogprices, vTime , vTesttime, iMeanBandwidth, iVarBandwidth, iPreAverage, iAcLag )
#endif
    
    return(lOut);
  
  
}
