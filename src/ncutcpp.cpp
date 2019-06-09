#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// [[Rcpp::export]]

double ncut_x(NumericVector x, double h, int n, int kmin){ // x must be sorted in increasing order
  NumericVector l(n);
  NumericVector r(n);
  l[0] = 1.0;// = r[n-1] = 1.0;
  for(int i=1; i<n; i++){
    l[i] = 1.0 + exp((x[i-1]-x[i])/h)*l[i-1];
    r[n-i-1] = exp((x[n-i-1]-x[n-i])/h)*(r[n-i]+1.0);
  }
  NumericVector cut(n);
  NumericVector vol(n);
  cut[0] = l[0]*r[0]/h;
  vol[0] = (l[0]+r[0])/h;//-1.0;
  for(int i=1; i<n; i++){
    cut[i] = l[i]*r[i]/h;
    vol[i] = vol[i-1] + (l[i]+r[i])/h;//-1.0;
  }
  double fval;
  double optval = 10.0;
  for(int i=kmin-1; i<(n-kmin); i++){
    fval = cut[i]*(1/(vol[i])+1/(vol[n-1]-vol[i]));
    if(fval < optval){
      optval = fval;
    }
  }
  //NumericVector output(n-1);
  //for(int i=0; i<(n-1); i++){
  //  output[i] = cut[i]*(1/(vol[i])+1/(vol[n-1]-vol[i]));
  //}
  //return output;
  return optval;
}



// [[Rcpp::export]]

NumericVector dncut_x(NumericVector x, double h, int n, int kmin){ // x must be sorted in increasing order
  NumericVector l(n);
  NumericVector r(n);
  l[0] = 1.0;// = r[n-1] = 1.0;
  for(int i=1; i<n; i++){
    l[i] = 1 + exp((x[i-1]-x[i])/h)*l[i-1];
    r[n-i-1] = exp((x[n-i-1]-x[n-i])/h)*(r[n-i]+1.0);
  }
  int k=kmin;
  double optval = 10.0;
  NumericVector cut(n);
  NumericVector vol(n);
  cut[0] = l[0]*r[0]/h;
  vol[0] = (r[0]+l[0])/h;//-1.0;
  for(int i=1; i<n; i++){
    cut[i] = l[i]*r[i]/h;
    vol[i] = vol[i-1] + (l[i]+r[i])/h;//-1.0;
  }
  double fval;
  for(int i=kmin-1; i<(n-kmin); i++){
    fval = cut[i]*(1/(vol[i])+1/(vol[n-1]-vol[i]));
    if(fval < optval){
      optval = fval;
      k = i+1;
    }
  }
  NumericVector dcut(n);
  NumericVector dvol(n);
  NumericVector dvoln(n);
  for(int i=0; i<k; i++){
    dcut[i] = exp((x[i]-x[k-1])/h)*r[k-1]/h/h;
    dvol[i] = (2.0*r[i]-exp((x[i]-x[k-1])/h)*r[k-1]-2.0*l[i])/h/h;
    dvoln[i] = 2.0/h/h*(r[i]-l[i]);
  }
  for(int i=k; i<n; i++){
    dcut[i] = -exp((x[k-1]-x[i])/h)*l[k-1]/h/h;
    dvol[i] = -exp((x[k-1]-x[i])/h)*l[k-1]/h/h;
    dvoln[i] = 2.0/h/h*(r[i]-l[i]);
  }
  NumericVector output(n);
  double ocut = cut[k-1];
  double ovol = vol[k-1];
  double onvol = vol[n-1]-vol[k-1];
  for(int i=0; i<n; i++){
    output[i] = dcut[i]*(1.0/ovol+1.0/onvol)-ocut*(dvol[i]/ovol/ovol+(dvoln[i]-dvol[i])/onvol/onvol);
  }
  return output;
}

