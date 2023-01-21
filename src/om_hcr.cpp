#include <TMB.hpp>

/*
 * Harvest strategies for fisheries with highly variable recruitment dynamics
 *                      Cahill & Walters oct 2022
 */

// linear 
template <class Type> 
Type ut_linear(vector<Type> par, Type vulb)
{ 
  vector<Type> seq(2); 
  seq.setZero(); 
  seq(0) = 0; 
  seq(1) = par(0)*(vulb-par(1))/vulb; 
  Type out = max(seq);  
  return out;
}  

// spline 
template <class Type> 
Type ut_spline(vector<Type> par, vector<Type> knots, Type vulb)
{ 
  tmbutils::splinefun<Type> sp(knots,par);
  Type out = sp(vulb);              
  return out;
}

// rectilinear 
template <class Type> 
Type ut_rect(vector<Type> par, Type vulb)
{ 
  vector<Type> seq1(2); 
  vector<Type> seq2(2); 
  seq1.setZero(); seq2.setZero();  
  
  seq1(0) = 0; 
  seq1(1) = par(0)*(vulb-par(1))/(par(2) - par(1)); 
  seq2(0) = par(0); 
  seq2(1) = max(seq1);
  
  Type out = min(seq2); 
  return out;
}

// exponential form 
template <class Type>
Type ut_exp(vector<Type> par, Type vulb)
{
  Type out = par(0)*(1-exp(-par(1)*vulb)); 
  return out;
}

// logit form 
template <class Type>
Type ut_logit(vector<Type> par, Type wbar, Type vulb)
{
  Type out = invlogit(par(0) + par(1)*vulb + par(2)*wbar); 
  return out;
}

// logit-linear
template <class Type>
Type ut_logit_linear(vector<Type> par, Type wbar, Type vulb)
{
  vector<Type> seq(2);
  seq.setZero();
  seq(0) = 0;
  seq(1) = par(3)*(vulb - par(0))/vulb/(1 + exp(-par(1)*(wbar - par(2)))); 
  Type out = max(seq);
  return out; 
}

// canada rule -- note dfopar read in as data
template <class Type> 
Type ut_dfo(vector<Type> dfopar, Type vulb)
{ 
  Type Umsy = dfopar(0);
  Type Bo = dfopar(1);
  Type blrp = 0.3*Bo; 
  Type ulrp = 0.5*Bo; 
  Type out = Umsy * (vulb - blrp) / (ulrp - blrp); 
  if (out < 0){out = 0;}
  if (out > Umsy){out = Umsy;}
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
  DATA_SCALAR(upow);        // utility power, 1 = risk neutral utility for catch, < 1 = risk aversion
  DATA_VECTOR(ages);    
  DATA_VECTOR(recmult);     // recruitment sequence
  DATA_INTEGER(hcrmode);    // feedback policy
  DATA_VECTOR(knots);       // spline knots
  DATA_VECTOR(dfopar);      // Umsy, Bmsy
  DATA_VECTOR(vmult);       // survey error = exp(sd_survey * (N(0,1)) - 0.5 * (sd_survey)^2)
  DATA_VECTOR(useq);        // sequence from 0 to 1 for stochastic Umsy estimation
  DATA_INTEGER(modulus);    // how often do we want to collapse the fishery?
  DATA_INTEGER(usequota);   // use a quota?
  DATA_INTEGER(umax);       // cap on implemented ut
  DATA_VECTOR(umult);       // implementation error ut
  DATA_SCALAR(dev);         // for smooth ut implementation
  
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

  for(int a = 0; a < n_age; a ++){
    vul(a) = 1 /( 1 + exp(-asl*(ages(a) - ahv))); 
    wt(a) = pow((1 - exp(-vbk*(ages(a)))), 3); 
    mat(a) = 1/(1 + exp(-asl*(ages(a) - ahm))); 
    if(a == 0){ 
      Lo(a) = 1;
      Lf(a) = 1; 
    }  
    else if(a > 0 && a < (n_age - 1)){
      Lo(a) = Lo(a-1)*s;
      Lf(a) = Lf(a-1)*s*(1 - vul(a-1)*uo);
    }
    else if(a == (n_age - 1)){          
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
  vector<Type> wbar(n_year);
  vector<Type> yield(n_year);
  vector<Type> utility(n_year);
  vector<Type> rec(n_year);
  vector<Type> ssb(n_year);
  vector<Type> vulb(n_year);
  vector<Type> vbobs(n_year);
  vector<Type> ut(n_year);
  vector<Type> tac(n_year);
  vector<Type> uout(n_year); 
  vector<Type> ftt(n_year); 
  vector<Type> ut2(n_year); 
  
  
  abar.setZero();wbar.setZero(); yield.setZero(); 
  utility.setZero(); rec.setZero(); ssb.setZero(); 
  vulb.setZero(); ut.setZero(); vbobs.setZero(); 
  tac.setZero(); uout.setZero(); ftt.setZero();
  ut2.setZero(); 
  
  n = ninit; 
  Type obj = 0;

  for(int t = 0; t < n_year; t++){
    if(t%modulus==0){n = rinit*n;}
    vulb(t) = (vul*n*wt).sum();  
    vbobs(t) = vulb(t)*vmult(t);  
    ssb(t) = (mwt*n).sum();                                          
    abar(t) = (vul*ages*n).sum() / sum(n);                             
    wbar(t) = (vul*n*wt).sum() / (n*wt).sum(); 
    switch(hcrmode){
      case 0:
        ut(t) = par(t);
      break;
        
      case 1:
        ut(t) = ut_linear(par, vbobs(t));
      break;
      
      case 2:
        ut(t) = ut_spline(par, knots, vbobs(t));
      break;

      case 3:
        ut(t) = ut_rect(par, vbobs(t));
      break;
      
      case 4:
        ut(t) = ut_exp(par, vbobs(t));
      break;
      
      case 5:
        ut(t) = ut_logit(par, wbar(t), vbobs(t));
      break;
      
      case 6:
        ut(t) = ut_logit_linear(par, wbar(t), vbobs(t));
      break;
      
      case 7:
        ut(t) = ut_dfo(dfopar, vbobs(t));
      break;
      
      default:
        std::cout<<"HCR code not yet implemented."<<std::endl;
      exit(EXIT_FAILURE);
      break;
    }
    if(usequota && hcrmode > 0){ // if usequota and hcrmode != OM
      tac(t) = ut(t)*vbobs(t);
      ut(t) = tac(t)/vulb(t);
      uout(t) = dev*log(exp((umax/dev)) + 1) - dev*log(exp(-(ut(t) - umax)/dev) + 1);
      ftt(t) = -log(1.0001 - uout(t))*(1 + umult(t));
      ut(t)= 1 - exp(-ftt(t));
    }
    yield(t) = ut2(t)*vulb(t);; 
    utility(t) = pow(yield(t), upow);
    n = s*n*(1-vul*ut(t)); 
    n(n_age - 1) = n(n_age - 1) + n(n_age - 2);                    
    for(int a = (n_age - 2); a > 0; a--){n(a) = n(a - 1);}        
    n(0) = reca*ssb(t)/(1 + recb*ssb(t))*recmult(t);   
    rec(t) = n(0); 
  }
  obj -= utility.sum();
  
  // stochastic estimation of umay, bo, bmay
  SIMULATE{ 
    Type umay = 0; 
    Type bmay = 0; 
    Type may = 0; 
    Type bo = 0; 
    for(int i = 0; i < useq.size(); i++){
      Type ut = useq(i); 
      n = ninit; 
      for(int t = 0; t < n_year; t++){
        vulb(t) = (vul*n*wt).sum();  
        ssb(t) = (mwt*n).sum();                                          
        yield(t) = ut*vulb(t);                                      
        n = s*n*(1-vul*ut); 
        n(n_age - 1) = n(n_age - 1) + n(n_age - 2);                    
        for(int a = (n_age - 2); a > 0; a--){n(a) = n(a - 1);}        
        n(0) = reca*ssb(t) / (1 + recb*ssb(t))*recmult(t); 
      }
      Type avg_yield = yield.sum()/yield.size();
      if(ut == 0){ // unfished state
        bo = vulb.sum()/vulb.size(); 
      } else if (avg_yield > may){  
        may = avg_yield; 
        umay = ut; 
        bmay = may / umay; 
      }
    }
    REPORT(may); 
    REPORT(umay); 
    REPORT(bmay); 
    REPORT(bo); 
  }
  
  REPORT(rec); 
  REPORT(ssb);
  REPORT(yield);
  REPORT(vulb);
  REPORT(vbobs); 
  REPORT(abar);
  REPORT(wbar); 
  REPORT(utility);
  REPORT(ut); 
  REPORT(tac);
  REPORT(uout);
  REPORT(ftt);
  REPORT(ut2); 

  return obj; 
}
