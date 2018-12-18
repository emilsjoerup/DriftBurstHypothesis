#include <RcppArmadillo.h>
#if defined(_OPENMP)
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace arma;

// [[Rcpp::export]]

arma::vec HACWeightC(int iLag){ 
  vec vW = linspace(1, iLag, iLag) / iLag;
  uvec vIdx = find(vW<=0.5);
  uvec vNotIdx = find(vW>0.5);
  vW.elem(vIdx) = 1-6*pow(vW.elem(vIdx),2) + 6 * pow(vW.elem(vIdx), 3);
  vW.elem(vNotIdx) = 2 * pow((1 - vW.elem(vNotIdx)), 3);
  return(vW);
  // Input will always be non-zero and positive, thus, the 0 case is ignored.
}

// [[Rcpp::export]]
int AutomaticLagSelectionC(arma::vec vX, double dMu){

  double dC = 2.6614;
  int iQ = 2;
  int iN = round(4.0 * pow((dMu/100), 0.16));
  
  double dRoot = 0.2;
  int iT = vX.size();
  vec vE = zeros(iN+1);
  vec vAc = zeros(iN);
  vec vFoo = linspace(1, (iN+1), (iN+1));
  
  for(int i=0; i<(iN+1); i++){
    vE[i] = mean(vX(span((i+1), iT-1)) % vX(span(0, (iT-i-2))) );
  }

  double d0 = - sum(vFoo % vE);
  
  for(int i=1; i<(iN+1); i++){
    vec vFoo = linspace(1, (iN-i+1), (iN-i+1));
    vAc[(i-1)] = sum(vFoo % vE(span((i), (iN))));
  }
  
  double s0 = d0 + iQ * sum( - vAc);
  
  
  double sQ = iQ * sum(pow(linspace(1, iN, iN), 2.0) % vAc);
  double dGamma = dC * pow(pow((sQ / s0), 2.0 ), dRoot);
  
  int iOut = round(dGamma * pow(dMu, dRoot));
  return(iOut);

  
}

//[[Rcpp::export]]
double AsymptoticVarianceC(arma::vec vIn, int iLag){
try{
  double d0 = sum(square(vIn));
  if(iLag == 0){
    vec vAc = zeros(1);
  }
  vec vAc = zeros(iLag);
  int iT = vIn.size();
  for(int i=0; i<iLag; i++){

    vAc[i] = sum(vIn(span((i+1), iT-1)) % vIn(span(0, (iT-i-2))));
  }
 vec vW = HACWeightC(iLag);
double dOut = d0 + 2.0 * sum(vAc % vW);
return(dOut);
}
catch( std::exception &ex ) {
  return(-1.0);
} catch(...) {
  return(-1.0);
}

}



//[[Rcpp::export]]

Rcpp::List DriftBurstLoopC(arma::vec vPreAveraged, arma::vec diffedlogprices, arma::vec vTime, arma::vec vTesttime, 
                          int iMeanBandwidth, int iVarBandwidth, int iPreAverage, int iAcLag){
  int iT = vTesttime.size();
  int iN = vPreAveraged.size();
  int iQ;
  double dReturnINF = 0.0;
  
  arma::vec vFooMu(iT);
  arma::vec vMu(iT);
  arma::vec vFooSigma(iT);
  arma::vec vSigma(iT);
  arma::vec vAcLag(iT);
  arma::uvec vIdx(iT);
  arma::vec vX = zeros<vec>(iN);
  arma::vec vWm = zeros<vec>(iN);
  arma::vec vWvar = zeros<vec>(iN);
  int iMaxLag = 20;
  for(int i=1; i<iT; i++){
    
    vX = vTime(span(0, (iN-1))) - vTesttime[i];
    vIdx = find(vX<=0);
    vWm = exp(-abs(vX.elem(vIdx)/iMeanBandwidth));  // exponential kernel
    vFooMu[i] = sum(vWm);
    
    
    vMu[i] = sum(vWm.elem(vIdx) % vPreAveraged.elem(vIdx)) / iMeanBandwidth; 
    
    vWvar = exp(-abs(vX.elem(vIdx)/iVarBandwidth));  // exponential kernel
  
    vFooSigma[i] = sum(vWvar);

    if(iAcLag == -1){
      iQ = AutomaticLagSelectionC(diffedlogprices, vFooMu[i]);
      vec vFoo(2); vFoo[0] = iQ ; vFoo[1] = iMaxLag;
      vAcLag[i] = min(vFoo) + 2 * (iPreAverage-1);
      
      
      double dAsympVar = AsymptoticVarianceC((vWvar.elem(vIdx) % vPreAveraged.elem(vIdx)),vAcLag[i]) / iVarBandwidth; 
      // currently this function sometimes(rarely) fails, thus a check is incorporated to prevent R from crashing
      
      if(dAsympVar == -1.0){
        dReturnINF = 1;
      }
      vSigma[i] = dAsympVar;
      
    }else{
      double dAsympVar = AsymptoticVarianceC((vWvar.elem(vIdx) % vPreAveraged.elem(vIdx)),iAcLag) / iVarBandwidth;
      vSigma[i] = dAsympVar;
    }


  }

  if(dReturnINF == 1){
    arma::vec vFoo(1);
    vFoo[0] = R_PosInf;
    return(Rcpp::List::create(Rcpp::Named("DriftBursts") = vFoo));
  }
  arma::vec vDb = sqrt(iMeanBandwidth) * vMu / sqrt(vSigma) ; 
  vDb[0] = 0;
  return(Rcpp::List::create(Rcpp::Named("DriftBursts") = vDb,
                            Rcpp::Named("Sigma") = vSigma,
                            Rcpp::Named("Mu") = vMu));
}



//[[Rcpp::export]]


Rcpp::List DriftBurstLoopCPAR(arma::vec vPreAveraged, arma::vec diffedlogprices, arma::vec vTime, 
                             arma::vec vTesttime, 
                             int iMeanBandwidth, int iVarBandwidth, int iPreAverage, int iAcLag , int iCores){
  
  
#if defined(_OPENMP)
  omp_set_num_threads(iCores);
  int iT = vTesttime.size();
  int iN = vPreAveraged.size();
  
  
  arma::vec vFooMu= zeros<vec>(iT);
  arma::vec vMu= zeros<vec>(iT);
  arma::vec vFooSigma= zeros<vec>(iT);
  arma::vec vSigma= zeros<vec>(iT);
  arma::uvec vIdx= zeros<uvec>(iT);
  arma::vec vX= zeros<vec>(iN);
  arma::vec vWm = zeros<vec>(iN);
  arma::vec vWvar = zeros<vec>(iN);
  
  double dReturnINF = 0.0;
  int i;
  int iMaxLag = 20;
  //Parallelization setup
  #pragma omp parallel for default(none)\
  shared(vPreAveraged, diffedlogprices, vTime, vTesttime, iMeanBandwidth, iVarBandwidth, iPreAverage, iAcLag, vMu, vFooMu, vFooSigma, vSigma,iT,iMaxLag, dReturnINF, iN) \
    private(vX, vWm, vWvar, i, vIdx)
    
    for(i=1; i<iT; i++){
      
      vX = vTime(span(0, (iN-1))) - vTesttime[i];
      vIdx = find(vX<=0);
      vWm = exp(-abs(vX.elem(vIdx)/iMeanBandwidth));  // exponential kernel
      vFooMu[i] = sum(vWm);
      
      vMu[i] = sum(vWm.elem(vIdx) % vPreAveraged.elem(vIdx)) / iMeanBandwidth; 
      vWvar = exp(-abs(vX.elem(vIdx)/iVarBandwidth));  // exponential kernel
      vFooSigma[i] = sum(vWvar);
      
      if(iAcLag == -1){
        int iQ = AutomaticLagSelectionC(diffedlogprices, vFooMu[i]);
        arma::vec vFoo(2); vFoo[0] = iQ ; vFoo[1] = iMaxLag;
        int iAutoAcLag = min(vFoo) + 2.0 * (iPreAverage-1);
        
        double dAsympVar = AsymptoticVarianceC((vWvar.elem(vIdx) % vPreAveraged.elem(vIdx)),iAutoAcLag) / iVarBandwidth;
        // currently this function rarely (0.5% of times (might be due to my data)) fails, 
        // thus a check is incorporated to prevent R from crashing
        if(dAsympVar == -1.0){
        dReturnINF = 1;
        }
        vSigma[i] = dAsympVar;
        
      }else{
        double dAsympVar = AsymptoticVarianceC((vWvar.elem(vIdx) % vPreAveraged.elem(vIdx)),iAcLag) / iVarBandwidth;
        vSigma[i] = dAsympVar;
      }
      
      
    }
  arma::vec vDb = sqrt(iMeanBandwidth) * vMu / sqrt(vSigma) ; 
#else
  Rf_warning("OpenMP is not available. Sequential processing is used.");
  vDb = DriftBurstLoopC(vPreAveraged, diffedlogprices, vTime , vTesttime, iMeanBandwidth, iVarBandwidth, iPreAverage, iAcLag )
#endif
  if(dReturnINF == 1){
    arma::vec vFoo(1);
    vFoo[0] = R_PosInf;
    return(Rcpp::List::create(Rcpp::Named("DriftBursts") = vFoo));
  }

  return(Rcpp::List::create(Rcpp::Named("DriftBursts") = vDb,
                            Rcpp::Named("Sigma") = vSigma,
                            Rcpp::Named("Mu") = vMu));

  
}
