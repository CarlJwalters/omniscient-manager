#include <TMB.hpp>
/*
 * omniscient manager for highly variable recruitment
 *             cahill & walters oct 2022
 */

template <class Type> 
Type ut_linear(vector<Type> par, Type vulb)
{ 
  // par(0) = ln(cslope), par(1) = ln(bmin)
  vector<Type> seq(2); 
  seq.setZero(); 
  seq(0) = 0; 
  seq(1) = exp(par(0))*(vulb-exp(par(1)))/vulb; 
  Type out = max(seq);               
  return out;
}  

template <class Type> 
Type ut_logistic(vector<Type> par, Type vulb)
{ 
  // Umax = par(0); bslope = par(1); bhalf = par(3); 
  Type out = par(0) / (1 + exp(-par(1)*(vulb - par(2)))); 
  return out;
}  

template <class Type> 
Type ut_spline(vector<Type> par, vector<Type> knots, Type vulb)
{ 
  tmbutils::splinefun<Type> sp(knots,par);
  Type TAC = sp(vulb); 
  Type out = TAC/vulb; 
  if(out < 0){ out = 0;} else if (out > 1){ out = 1;}
  return out;
}

template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_INTEGER(n_year); 
  DATA_INTEGER(n_age);
  DATA_SCALAR(vbk);         // von bertalanffy k 
  DATA_SCALAR(s);           // avg survival
  DATA_SCALAR(cr);          // compensation ratio
  DATA_SCALAR(rinit);       // initial number age 1 recruits
  DATA_SCALAR(ro);          // average unfished recruitment
  DATA_SCALAR(uo);          // average historical exploitation rate
  DATA_SCALAR(asl);         // vul parameter 1
  DATA_SCALAR(ahv);         // vul parameter 2
  DATA_SCALAR(ahm);         // age half mature 
  DATA_SCALAR(upow);        // utility power
  DATA_VECTOR(ages);    
  DATA_VECTOR(recmult);     // recruitment sequence
  DATA_INTEGER(objmode);    // 0 = MAY, 1 = HARA utility
  DATA_INTEGER(hcrmode);    // rule 0 = U(t); 1 = linear; 2 = logistic; 3 = spline  
  DATA_VECTOR(knots); 

  vector<Type> n(n_age);
  vector<Type> ninit(n_age);
  vector<Type> vul(n_age);
  vector<Type> wt(n_age);
  vector<Type> mat(n_age);
  vector<Type> Lo(n_age);   
  vector<Type> mwt(n_age);
  vector<Type> Lf(n_age);
  n.setZero(); ninit.setZero(); vul.setZero(); wt.setZero(); mat.setZero();
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
    else if(a > 0 && a < (n_age - 1)){ // ages 2-19
      Lo(a) = Lo(a-1)*s;
      Lf(a) = Lf(a-1)*s*(1 - vul(a-1)*uo);
    }
    else if(a == (n_age - 1)){         // plus group age 20
      Lo(a) = Lo(a - 1)*s / (1 - s); 
      Lf(a) = Lf(a - 1)*s*(1 - vul(a-1)*uo) / (1 - s*(1 - vul(a-1)*uo)); 
    }
  } 
  
  ninit = rinit*Lf; 
  mwt = mat*wt; 
  sbro = (Lo*mwt).sum(); 
  Type reca = cr/sbro; 
  Type recb = (cr - 1) / (ro*sbro); 
  
  PARAMETER_VECTOR(par); 
  
  vector<Type> abar(n_year);
  vector<Type> yield(n_year);
  vector<Type> utility(n_year);
  vector<Type> ssb(n_year);
  vector<Type> vulb(n_year);
  vector<Type> ut(n_year);
  abar.setZero(); yield.setZero(); utility.setZero();
  ssb.setZero(); vulb.setZero(); ut.setZero(); 
  
  n = ninit; 
  Type obj = 0;
  
  for(int t = 0; t < n_year; t++){
    if(t%100==0){n = ninit;}
    vulb(t) = (vul*n*wt).sum();                                    
    ssb(t) = (mwt*n).sum();                                        
    abar(t) = (ages*n).sum() / sum(n);                             
    
    switch(hcrmode){
      case 0:
        ut(t) = par(t);
      break;
        
      case 1:
        ut(t) = ut_linear(par, vulb(t));
      break;
        
      case 2:
        ut(t) = ut_logistic(par, vulb(t));
      break;
      
      case 3:
        ut(t) = ut_spline(par, knots, vulb(t));
      break;

      default:
        std::cout<<"Ut code not yet implemented."<<std::endl;
      exit(EXIT_FAILURE);
      break;
    }
    
    yield(t) = ut(t)*vulb(t);                                      
    utility(t) = pow(yield(t), upow);
    n = s*n*(1-vul*ut(t)); 
    n(n_age - 1) = n(n_age - 1) + n(n_age - 2);                    
    for(int a = (n_age - 2); a > 0; a--){n(a) = n(a - 1);}        
    n(0) = reca*ssb(t) / (1 + recb*ssb(t))*recmult(t);             
    
    switch(objmode){
      case 0:
        obj -= yield(t);
      break;
      
      case 1:
        obj -= utility(t);
      break;
      
      default:
        std::cout<<"Objective code not yet implemented."<<std::endl;
      exit(EXIT_FAILURE);
      break;
    }
  }
  
  REPORT(ssb);
  REPORT(yield);
  REPORT(vulb);
  REPORT(abar);
  REPORT(utility);
  REPORT(ut); 

  return obj; 
}
