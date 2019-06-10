#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


double d_abs(double x){
  double ret;
  if(x>0) ret = x;
  else ret = -x;
  return ret;
}

// [[Rcpp::export]]

double f_md_cpp(arma::vec v, arma::mat X, int n, int d, double h, double al, double C){
  double nv2=0;
  for(int i=0; i<d; i++) nv2 += v[i]*v[i];
  double nv = pow(nv2, 0.5);
  arma::vec v_(d);
  for(int i=0; i<d; i++) v_[i] = v[i]/nv;
  arma::vec p = X*v_;
  double miny = 0.0;
  if(al>0.00001){
    double var = 0;
    for(int i=0; i<n; i++) var += p[i]*p[i];
    var/=(n-1.0);
    double sd = pow(var, .5);
    NumericVector op(n);
    for(int i=0; i<n; i++) op[i] = p[i];
    std::sort(op.begin(),op.end());

    NumericMatrix L(2, n);
    NumericMatrix R(2, n);
    for(int i=0; i<=1; i++) L(i,0) = pow(-op[0], i);
    for(int i=1; i<n; i++){
      for(int j=0; j<=1; j++){
        L(j,i) = pow(-op[i],j) + exp((op[i-1]-op[i])/h)*L(j,i-1);
        R(j,n-i-1) = exp((op[n-i-1]-op[n-i])/h)*(pow(op[n-i],j)+R(j,n-i));
      }
    }
    NumericVector df(n);
    double denom = 4.0*n*pow(h, 3);
    for(int i=0; i<n; i++){
      for(int j=0; j<=1; j++) df[i] -= (pow(op[i], 1-j)*L(j,i)-pow(-op[i],1-j)*R(j,i))/denom;
      if(op[i]<(-al*sd)) df[i] -= 2.0*C*(-al*sd-op[i]);
      if(op[i]>(al*sd)) df[i] += 2.0*C*(op[i]-al*sd);
    }
    double f_at_min, df_at_min;
    miny = 100000;
    double minx;
    int pos = 0;
    double lo, hi, mid;
    double eps = pow(0.1, 8);
    double exp_mult;
    while(pos<(n-1)){
      if(df[pos]<0 && df[pos+1]>0){
        denom = 4.0*n*pow(h, 3);
        f_at_min = 0.0;
        df_at_min = 1.0;
        lo = op[pos];
        hi = op[pos+1];
        mid = 0.5*lo+0.5*hi;
        while((hi-lo)>eps && d_abs(df_at_min)>eps){
          df_at_min = 0.0;
          mid = 0.5*lo+0.5*hi;
          exp_mult = exp((op[pos]-mid)/h);
          for(int j=0; j<=1; j++) df_at_min -= (pow(mid,1-j)*exp_mult*L(j,pos)-pow(-mid,1-j)/exp_mult*R(j,pos))/denom;
          if(mid<(-al*sd)) df_at_min -= 2.0*C*(-al*sd-mid);
          if(mid>(al*sd)) df_at_min += 2.0*C*(mid-al*sd);
          if(df_at_min<(-eps)) lo = mid;
          if(df_at_min>eps) hi = mid;
        }
        exp_mult = exp((op[pos]-mid)/h);
        for(int orddo = 0; orddo<=1; orddo++){
          denom = 4.0*n*pow(h,orddo+1);
          for(int j=0; j<=orddo; j++) f_at_min += (pow(mid, orddo-j)*L(j,pos)*exp_mult+pow(-mid,orddo-j)*R(j,pos)/exp_mult)/denom;
        }
        if(mid<(-al*sd)) f_at_min += C*pow(-al*sd-mid,2);
        if(mid>(al*sd)) f_at_min += C*pow(mid-al*sd,2);
        if(f_at_min < miny){
          miny = f_at_min;
          minx = mid;
        }
      }
      pos += 1;
    }
  }
  else{
    double constnt = 1.0/n/h/4.0;
    for(int i=0; i<n; i++) miny += exp(-d_abs(p[i])/h)*(1.0 + d_abs(p[i])/h)*constnt;
  }
  return(miny);
}

// [[Rcpp::export]]

arma::vec df_md_cpp(arma::vec v, arma::mat X, int n, int d, double h, double al, double C){
  double nv2=0;
  for(int i=0; i<d; i++) nv2 += v[i]*v[i];
  double nv = pow(nv2, 0.5);
  arma::vec v_(d);
  for(int i=0; i<d; i++) v_[i] = v[i]/nv;
  arma::vec p = X*v_;
  double var = 0;
  for(int i=0; i<n; i++) var += p[i]*p[i];
  var/=(n-1.0);
  double sd = pow(var, .5);
  double minx = 0;
  double denom;
  if(al > 0.00001){
    NumericVector op(n);
    for(int i=0; i<n; i++) op[i] = p[i];
    std::sort(op.begin(),op.end());

    NumericMatrix L(2, n);
    NumericMatrix R(2, n);
    for(int i=0; i<=1; i++) L(i,0) = pow(-op[0], i);
    for(int i=1; i<n; i++){
      for(int j=0; j<=1; j++){
        L(j,i) = pow(-op[i],j) + exp((op[i-1]-op[i])/h)*L(j,i-1);
        R(j,n-i-1) = exp((op[n-i-1]-op[n-i])/h)*(pow(op[n-i],j)+R(j,n-i));
      }
    }
    NumericVector df(n);
    denom = 4.0*n*pow(h, 3);
    for(int i=0; i<n; i++){
      for(int j=0; j<=1; j++) df[i] -= (pow(op[i], 1-j)*L(j,i)-pow(-op[i],1-j)*R(j,i))/denom;
      if(op[i]<(-al*sd)) df[i] -= 2.0*C*(-al*sd-op[i]);
      if(op[i]>(al*sd)) df[i] += 2.0*C*(op[i]-al*sd);
    }
    double f_at_min, df_at_min;
    double miny = 100000;
    int pos = 0;
    double lo, hi, mid;
    double eps = pow(0.1, 8);
    double exp_mult;
    while(pos<(n-1)){
      if(df[pos]<0 && df[pos+1]>0){
        denom = 4.0*n*pow(h, 3);
        lo = op[pos];
        hi = op[pos+1];
        mid = 0.5*lo+0.5*hi;
        df_at_min = 1.0;
        f_at_min = 0;
        while((hi-lo)>eps && d_abs(df_at_min)>eps){
          df_at_min = 0;
          mid = 0.5*lo+0.5*hi;
          exp_mult = exp((op[pos]-mid)/h);
          for(int j=0; j<=1; j++) df_at_min -= (pow(mid,1-j)*exp_mult*L(j,pos)-pow(-mid,1-j)/exp_mult*R(j,pos))/denom;
          if(mid<(-al*sd)) df_at_min -= 2.0*C*(-al*sd-mid);
          if(mid>(al*sd)) df_at_min += 2.0*C*(mid-al*sd);
          if(df_at_min<(-eps)) lo = mid;
          if(df_at_min>eps) hi = mid;
        }
        exp_mult = exp((op[pos]-mid)/h);
        for(int orddo = 0; orddo<=1; orddo++){
          denom = 4.0*n*pow(h,orddo+1);
          for(int j=0; j<=orddo; j++) f_at_min += (pow(mid, orddo-j)*L(j,pos)*exp_mult+pow(-mid,orddo-j)*R(j,pos)/exp_mult)/denom;
        }
        if(mid<(-al*sd)) f_at_min += C*pow(-al*sd-mid,2);
        if(mid>(al*sd)) f_at_min += C*pow(mid-al*sd,2);
        if(f_at_min < miny){
          miny = f_at_min;
          minx = mid;
        }
      }
      pos += 1;
    }
  }
  arma::vec dp(n);
  denom = 4.0*n*pow(h,3)*nv;
  for(int i=0; i<n; i++) dp[i] = (minx-p[i])*exp(-d_abs(minx-p[i])/h)/denom;
  if(minx<(-al*sd)){
    double cnst = 2.0*al*C/sd/(n-1.0)*(minx+al*sd)/nv;
    for(int i=0; i<n; i++) dp[i] += cnst*p[i];
  }
  if(minx>(al*sd)){
    double cnst = 2.0*al*C/sd/(n-1.0)*(al*sd-minx)/nv;
    for(int i=0; i<n; i++) dp[i] += cnst*p[i];
  }
  //arma::vec ret = dp.t()*X-(dp.t()*p)*v.t();
  double dpp = 0;
  for(int i=0; i<n; i++) dpp += dp[i]*p[i];
  arma::vec ret = X.t()*dp-dpp*v_;
  //arma::vec ret = X.t()*dp - v_*v_.t()*X.t()*dp;
  return(ret);
}


double norm_vec(arma::vec v, int d){
  double nm2 = 0;
  for(int i=0; i<d; i++) nm2 += v[i]*v[i];
  double ret = pow(nm2, 0.5);
  return(ret);
}


// [[Rcpp::export]]

int ismin_cpp(arma::vec v, arma::mat X, int n, int d, double h, double al, double C){
  double nv2=0;
  for(int i=0; i<d; i++) nv2 += v[i]*v[i];
  double nv = pow(nv2, 0.5);
  arma::vec v_(d);
  for(int i=0; i<d; i++) v_[i] = v[i]/nv;
  arma::vec p = X*v_;
  double var = 0;
  for(int i=0; i<n; i++) var += p[i]*p[i];
  var/=(n-1.0);
  double sd = pow(var, .5);
  NumericVector op(n);
  for(int i=0; i<n; i++) op[i] = p[i];
  std::sort(op.begin(),op.end());

  NumericMatrix L(2, n);
  NumericMatrix R(2, n);
  for(int i=0; i<=1; i++) L(i,0) = pow(-op[0], i);
  for(int i=1; i<n; i++){
    for(int j=0; j<=1; j++){
      L(j,i) = pow(-op[i],j) + exp((op[i-1]-op[i])/h)*L(j,i-1);
      R(j,n-i-1) = exp((op[n-i-1]-op[n-i])/h)*(pow(op[n-i],j)+R(j,n-i));
    }
  }
  NumericVector df(n);
  NumericVector df2(n);
  double denom = 4.0*n*pow(h, 3);
  for(int i=0; i<n; i++){
    for(int j=0; j<=1; j++){
      df[i] -= (pow(op[i], 1-j)*L(j,i)-pow(-op[i],1-j)*R(j,i))/denom;
      df2[i] -= (pow(op[i], 1-j)*L(j,i)-pow(-op[i],1-j)*R(j,i))/denom;
    }
    if(op[i]<(-al*sd)) df[i] -= 2*C*(-al*sd-op[i]);
    if(op[i]>(al*sd)) df[i] += 2*C*(op[i]-al*sd);
  }
  double f_at_min, df_at_min;
  double miny = 100000;
  double minx = 0;
  int pos = 0;
  double lo, hi, mid;
  double eps = pow(0.1, 7);
  double exp_mult;
  double mode1 = 0;
  double modef = 0;
  int mode1_found = 0;
  while(pos<(n-1)){
    if(df2[pos]>=0 && df2[pos+1]<=0){
      if(mode1_found==0){
        mode1_found = 1;
        mode1 = op[pos];
      }
      modef = op[pos];
    }
    if(df[pos]<=0 && df[pos+1]>=0){
      denom = 4.0*n*pow(h, 3);
      f_at_min = 0.0;
      df_at_min = 1.0;
      lo = op[pos];
      hi = op[pos+1];
      mid = 0.5*lo+0.5*hi;
      while((hi-lo)>eps && d_abs(df_at_min)>eps){
        df_at_min = 0.0;
        mid = 0.5*lo+0.5*hi;
        exp_mult = exp((op[pos]-mid)/h);
        for(int j=0; j<=1; j++) df_at_min -= (pow(mid,1-j)*exp_mult*L(j,pos)-pow(-mid,1-j)/exp_mult*R(j,pos))/denom;
        if(mid<(-al*sd)) df_at_min -= 2*C*(-al*sd-mid);
        if(mid>(al*sd)) df_at_min += 2*C*(mid-al*sd);
        if(df_at_min<(-eps)) lo = mid;
        if(df_at_min>eps) hi = mid;
      }
      exp_mult = exp((op[pos]-mid)/h);
      for(int orddo = 0; orddo<=1; orddo++){
        denom = 4.0*n*pow(h,orddo+1);
        for(int j=0; j<=orddo; j++) f_at_min += (pow(mid, orddo-j)*L(j,pos)*exp_mult+pow(-mid,orddo-j)*R(j,pos)/exp_mult)/denom;
      }
      if(mid<(-al*sd)) f_at_min += C*pow(-al*sd-mid,2);
      if(mid>(al*sd)) f_at_min += C*pow(mid-al*sd,2);
      if(f_at_min < miny){
        miny = f_at_min;
        minx = mid;
      }
    }
    pos += 1;
  }
  int ret = 0;
  if(minx>mode1 && minx<modef) ret = 1;
  return(ret);
}

// [[Rcpp::export]]

double md_b_cpp(arma::vec v, arma::mat X, int n, int d, double h, double al, double C){
  double nv2=0;
  for(int i=0; i<d; i++) nv2 += v[i]*v[i];
  double nv = pow(nv2, 0.5);
  arma::vec v_(d);
  for(int i=0; i<d; i++) v_[i] = v[i]/nv;
  arma::vec p = X*v_;
  double var = 0;
  for(int i=0; i<n; i++) var += p[i]*p[i];
  var/=(n-1.0);
  double sd = pow(var, .5);
  NumericVector op(n);
  for(int i=0; i<n; i++) op[i] = p[i];
  std::sort(op.begin(),op.end());

  NumericMatrix L(2, n);
  NumericMatrix R(2, n);
  for(int i=0; i<=1; i++) L(i,0) = pow(-op[0], i);
  for(int i=1; i<n; i++){
    for(int j=0; j<=1; j++){
      L(j,i) = pow(-op[i],j) + exp((op[i-1]-op[i])/h)*L(j,i-1);
      R(j,n-i-1) = exp((op[n-i-1]-op[n-i])/h)*(pow(op[n-i],j)+R(j,n-i));
    }
  }
  NumericVector df(n);
  double denom = 4.0*n*pow(h, 3);
  for(int i=0; i<n; i++){
    for(int j=0; j<=1; j++) df[i] -= (pow(op[i], 1-j)*L(j,i)-pow(-op[i],1-j)*R(j,i))/denom;
    if(op[i]<(-al*sd)) df[i] -= 2*C*(-al*sd-op[i]);
    if(op[i]>(al*sd)) df[i] += 2*C*(op[i]-al*sd);
  }
  double f_at_min, df_at_min;
  double miny = 100000;
  double minx = 0;
  int pos = 0;
  double lo, hi, mid;
  double eps = pow(0.1, 7);
  double exp_mult;
  while(pos<(n-1)){
    if(df[pos]<0 && df[pos+1]>0){
      denom = 4.0*n*pow(h, 3);
      f_at_min = 0.0;
      df_at_min = 1.0;
      lo = op[pos];
      hi = op[pos+1];
      mid = 0.5*lo+0.5*hi;
      while((hi-lo)>eps && d_abs(df_at_min)>eps){
        df_at_min = 0.0;
        mid = 0.5*lo+0.5*hi;
        exp_mult = exp((op[pos]-mid)/h);
        for(int j=0; j<=1; j++) df_at_min -= (pow(mid,1-j)*exp_mult*L(j,pos)-pow(-mid,1-j)/exp_mult*R(j,pos))/denom;
        if(mid<(-al*sd)) df_at_min -= 2*C*(-al*sd-mid);
        if(mid>(al*sd)) df_at_min += 2*C*(mid-al*sd);
        if(df_at_min<0) lo = mid;
        else hi = mid;
      }
      exp_mult = exp((op[pos]-mid)/h);
      for(int orddo = 0; orddo<=1; orddo++){
        denom = 4.0*n*pow(h,orddo+1);
        for(int j=0; j<=orddo; j++) f_at_min += (pow(mid, orddo-j)*L(j,pos)*exp_mult+pow(-mid,orddo-j)*R(j,pos)/exp_mult)/denom;
      }
      if(mid<(-al*sd)) f_at_min += C*pow(-al*sd-mid,2);
      if(mid>(al*sd)) f_at_min += C*pow(mid-al*sd,2);
      if(f_at_min < miny){
        miny = f_at_min;
        minx = mid;
      }
    }
    pos += 1;
  }
  return(minx);
}
