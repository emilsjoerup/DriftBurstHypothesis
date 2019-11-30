#include <RcppArmadillo.h>
#include "utils.h"
#if defined(_OPENMP)
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif
using namespace arma;
using namespace Rcpp;
const int iMaxLag = 20;
const double dC = 2.6614;


// [[Rcpp::export]]
arma::vec HACWeightC(int iLag){ 
  // Input will always be non-zero and positive, thus, we can ignore the 0 case and save some time
  vec vW = linspace(1, iLag,iLag) / iLag;
  int iIdx = floor(iLag/2.0);
  vW(span(0,iIdx-1)) = 1.0 - 6.0 * square(vW(span(0,iIdx-1))) + 6.0 * pow(vW(span(0,iIdx-1)), 3.0);
  vW(span(iIdx, iLag-1)) = 2.0 * pow(1 - vW(span(iIdx, iLag-1)), 3.0);
  return(vW);
}



//[[Rcpp::export]]

double AsymptoticVarianceC(arma::vec vIn, int iLag){  
  int iT = vIn.size();
  arma::vec vW;
  if(iT<=iLag){ //We cannot calculate the variance in this case as there are not enough observations.
    return(datum::nan);
  }else{
    vec vAc = zeros(iLag);
    if(iLag == 0 || iLag == 1){
      return(sum(vIn % vIn)); //Special cases where the weight will be 0.
    }else{
      vW = HACWeightC(iLag);
      for(int i=0; i<iLag; i++){
        vAc[i] = sum(vIn(span((i+1), iT-1)) % vIn(span(0, (iT-i-2))));
      }
    }
    double dOut = sum(vIn % vIn) + 2.0 * sum(vAc % vW);
    return(dOut);
  }
}



// [[Rcpp::export]]

double AutomaticLagSelectionC(arma::vec vX, double dMu){
  int iN = round(4.0 * pow((dMu/100), 0.16));
  int iT = vX.size();
  if(iT<=iN+1){
    return(iMaxLag);
  }
  
  vec vE = zeros(iN+1);
  vec vAc = zeros(iN);
  
  for(int i=0; i<=iN; i++){
    vE[i] = mean(vX(span((i+1), iT-1)) % vX(span(0, (iT-i-2))));
  }
  double d0 = - sum(linspace(1, (iN+1), (iN+1)) % vE);
  
  for(int i=1; i<=iN; i++){
    vAc[(i-1)] = sum(linspace(1, (iN-i+1), (iN-i+1)) % vE(span(i, iN)));
  }
  
  double s0 = d0 + 2.0 * -sum(vAc);
  
  
  double sQ = 2.0 * sum(square(linspace(1, iN, iN)) % vAc);
  double dGamma = dC * pow(pow((sQ / s0), 2.0 ), 0.2);
  
  int iOut = round(dGamma * pow(dMu, 0.2));
  
  return(iOut);
  
  
}


//[[Rcpp::export]]

Rcpp::List DriftBurstLoopC(arma::vec vPreAveraged, arma::vec diffedlogprices, arma::vec vTime, arma::vec vTesttime, 
                           double iMeanBandwidth, double iVarBandwidth, int iPreAverage, int iAcLag){
  int iT = vTesttime.size();
  int iQ;
  int iIdx;
  double iAutoAcLag;
  arma::vec vMu = zeros<vec>(iT);
  arma::vec vSigma = zeros<vec>(iT);
  arma::vec vX = zeros<vec>(vTime.size());
  arma::vec vWm;
  arma::vec vWvar;
  arma::vec vFoo(2);
  vFoo[1] = iMaxLag;
  
  for(int i = 0; i<iT; i++){
    vX = vTime - vTesttime[i];
    iIdx = sum((vX<=0))-2;
    
    vWm = exp(vX(span(0,iIdx)) / iMeanBandwidth);  // exponential kernel
    
    vMu[i] = sum(vWm(span(0,iIdx)) % vPreAveraged(span(0,iIdx))); 
    
    vWvar = exp(vX(span(0,iIdx)) / iVarBandwidth);  // exponential kernel
    
    if(iAcLag == -1){
      iQ = AutomaticLagSelectionC(diffedlogprices(span(0,iIdx)), sum(vWm));
      vFoo[0] = iQ;
      iAutoAcLag = min(vFoo) + 2.0 * (iPreAverage-1);
      vSigma[i] = AsymptoticVarianceC((vWvar(span(0,iIdx)) % vPreAveraged(span(0,iIdx))), iAutoAcLag); 
      
      
    }else{
      vSigma[i] = AsymptoticVarianceC((vWvar(span(0,iIdx)) % vPreAveraged(span(0,iIdx))), iAcLag);
  
    }
    
    
  }
  
  arma::vec vDB = sqrt(iMeanBandwidth) * (vMu / iMeanBandwidth) / sqrt(vSigma/ iVarBandwidth) ; 
  Rcpp::List lOut = Rcpp::List::create(Rcpp::Named("driftBursts") = vDB,
                                        Rcpp::Named("sigma") = vSigma / iVarBandwidth,
                                        Rcpp::Named("mu") = vMu / iMeanBandwidth);
  return(lOut);
}




//[[Rcpp::export]]

Rcpp::List DriftBurstLoopCPAR(arma::vec vPreAveraged, arma::vec diffedlogprices, arma::vec vTime, 
                              arma::vec vTesttime, double iMeanBandwidth, double iVarBandwidth, 
                              int iPreAverage, int iAcLag , int iCores){
#if defined(_OPENMP)
  omp_set_num_threads(iCores);
  int iT = vTesttime.size();
  int iIdx, iAutoAcLag;
  double iQ;
  arma::vec vMu= zeros<vec>(iT);
  arma::vec vSigma= zeros<vec>(iT);
  arma::vec vX = zeros<vec>(vTime.size());
  arma::vec vWvar = zeros<vec>(vTime.size());
  arma::vec vWm =  zeros<vec>(vTime.size());
  arma::vec vFoo = {0, iMaxLag};
  //Parallelization setup
  #pragma omp parallel for default(none)\
  shared(vPreAveraged, diffedlogprices, vTime, vTesttime, iMeanBandwidth,\
         iVarBandwidth, iPreAverage, iAcLag, vMu, vSigma, iT)\
    private(vX, vWm, vWvar, iIdx, iQ, iAutoAcLag, vFoo)
    
    for(int i = 0 ; i<iT; i++){
      
      vX = vTime - vTesttime[i];
      
      iIdx = sum((vX <= 0)) - 2;
      
      vWm = exp(vX(span(0, iIdx)) / iMeanBandwidth);  // exponential kernel
      
      vMu[i] = sum(vWm(span(0,iIdx)) % vPreAveraged(span(0, iIdx)));
      
      vWvar = exp(vX(span(0,iIdx)) / iVarBandwidth);  // exponential kernel
      
      if(iAcLag == -1){
        
        iQ = AutomaticLagSelectionC(diffedlogprices(span(0,iIdx)), sum(vWm));
        
        arma::vec vFoo = {iQ, iMaxLag};
        iAutoAcLag = min(vFoo) + 2.0 * (iPreAverage-1);
        vSigma[i]= AsymptoticVarianceC((vWvar(span(0,iIdx)) % vPreAveraged(span(0,iIdx))), iAutoAcLag); 
      }else{
        vSigma[i] = AsymptoticVarianceC((vWvar(span(0,iIdx)) % vPreAveraged(span(0,iIdx))), iAcLag);
      }
      
      
      
    }
    
  arma::vec vDB = sqrt(iMeanBandwidth) * (vMu/ iMeanBandwidth) / sqrt(vSigma/ iVarBandwidth) ; 
  Rcpp::List lOut =  Rcpp::List::create(Rcpp::Named("driftBursts") = vDB,
                                        Rcpp::Named("sigma") = vSigma / iVarBandwidth,
                                        Rcpp::Named("mu") = vMu / iMeanBandwidth);
#else
  // If openMP is not available, we use single core execution.
  Rf_warning("OpenMP is not available. Sequential processing is used.");
  Rcpp::List lOut = DriftBurstLoopC(vPreAveraged, diffedlogprices, vTime , vTesttime, iMeanBandwidth, iVarBandwidth, iPreAverage, iAcLag )
#endif
    
    return(lOut);
  
  
}

