#include <TMB.hpp>
/*
 * omniscient manager for highly variable recruitment
 *             cahill & walters oct 2022
 */
template <class Type> 
Type ut_logistic(vector<Type> par, Type vulb, Type xinc)
{ 
  int ix = CppAD::Integer(vulb / xinc);
  if(ix > par.size() - 2){ ix = par.size() - 2; }
  Type y1 = par(ix);
  Type y2 = par(ix + 1);
  Type out = y1 + (vulb/xinc - (ix))* (y2 - y1);
  return out;
}  

template <class Type> 
Type ut_linear(vector<Type> par, Type vulb)
{ 
  Type G = 10; 
  Type offset = 1e-6; 
  Type TAC = offset + (par(1)*(vulb - par(0)))/(1+exp(-G*(vulb-par(0)))); 
  Type out = TAC/(vulb + offset);                                           
  return out;
}  

template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_INTEGER(n_year); 
  DATA_INTEGER(n_age);
  DATA_SCALAR(vbk);      // von bertalanffy k 
  DATA_SCALAR(s);        // avg survival
  DATA_SCALAR(cr);       // compensation ratio
  DATA_SCALAR(rinit);    // initial number age 1 recruits
  DATA_SCALAR(ro);       // average unfished recruitment
  DATA_SCALAR(uo);       // average historical exploitation rate
  DATA_SCALAR(asl);      // vul parameter 1
  DATA_SCALAR(ahv);      // vul parameter 2
  DATA_SCALAR(ahm);      // age half mature 
  DATA_SCALAR(upow);     // utility power
  DATA_VECTOR(ages); 
  DATA_VECTOR(recmult);  // recruitment sequence
  DATA_INTEGER(obj_ctl); // 0 = MAY, 1 = HARA utility
  DATA_INTEGER(hcr);     // which harvest control rule 0 = linear; 1 = misery 
  
  vector<Type> n(n_age);
  vector<Type> vul(n_age);
  vector<Type> wt(n_age);
  vector<Type> mat(n_age);
  vector<Type> Lo(n_age);   
  vector<Type> mwt(n_age);
  vector<Type> Lf(n_age);
  n.setZero(); vul.setZero(); wt.setZero(); mat.setZero();
  Lo.setZero(); mwt.setZero(); Lf.setZero(); Type sbro = 0;  
  
  // set up leading vectors
  for(int a = 0; a < n_age; a ++){
    vul(a) = 1 /( 1 + exp(-asl*(ages(a) - ahv))); 
    wt(a) = pow((1 - exp(-vbk*(ages(a)))), 3); 
    mat(a) = 1/(1 + exp(-asl*(ages(a) - ahm))); 
    if(a == 0){ 
      Lo(a) = 1;
      Lf(a) = 1; 
    } 
    if(a > 0 && a < (n_age - 1)){ // ages 2-19
      Lo(a) = Lo(a-1)*s;
      Lf(a) = Lf(a-1)*s*(1 - vul(a-1)*uo);
    }
    if(a == (n_age - 1)){         // plus group age 20
      Lo(a) = Lo(a - 1)*s / (1 - s); 
      Lf(a) = Lf(a - 1)*s*(1 - vul(a-1)*uo) / (1 - s*(1 - vul(a-1)*uo)); 
    }
  } 
  n = rinit*Lf; 
  mwt = mat*wt; 
  sbro = (Lo*mwt).sum(); 
  Type reca = cr/sbro; 
  Type recb = (cr - 1) / (ro*sbro); 
  
  // parameters to solve 
  PARAMETER_VECTOR(par); 
  
  vector<Type> abar(n_year);
  vector<Type> yield(n_year);
  vector<Type> utility(n_year);
  vector<Type> ssb(n_year);
  vector<Type> vulb(n_year);
  vector<Type> ut(n_year);
  abar.setZero(); yield.setZero(); utility.setZero();
  ssb.setZero(); vulb.setZero(); ut.setZero(); 
  
  // run simulation 
  for(int t = 0; t < n_year; t++){
    vulb(t) = (vul*n*wt).sum();                                    // sumproduct(vul*n*w) across a
    ssb(t) = (mwt*n).sum();                                        // sumproduct(mwt * n)
    abar(t) = (ages*n).sum() / sum(n);                             // sumproduct(ages*n) / sum(n)
    if(hcr == 0){ut(t) = par(t);}
    if(hcr == 1){ut(t) = ut_linear(par, vulb(t));}                 // linear hcr
    //if(hcr == 1){ut(t) = ut_map(par, vulb(t), xinc);}              
    yield(t) = ut(t)*vulb(t);                                      
    utility(t) = pow(yield(t), upow);
    n = s*n*(1-vul*ut(t)); 
    n(n_age - 1) = n(n_age - 1) + n(n_age - 2);                    // plus group
    for(int a = (n_age - 2); a > 0; a--){n(a) = n(a - 1);}         // advance fish one a
    n(0) = reca*ssb(t) / (1 + recb*ssb(t))*recmult(t);             // recruits
  }
  
  REPORT(ssb);
  REPORT(yield);
  REPORT(vulb);
  REPORT(abar);
  REPORT(utility);
  REPORT(ut); 

  // objective function
  Type obj = 0;
  if(obj_ctl == 0){  // yield objective
    obj -= sum(yield);
  }
  if(obj_ctl == 1){  // hara utility objective
    obj -= sum(utility);
  }
  // if(hcr == 1){
  //  for(int i = 0; i < par.size(); i++){ obj += 0.000001*pow(par(i),2); }
  // }
  return obj; 
}

