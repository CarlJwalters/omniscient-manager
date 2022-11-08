#include <TMB.hpp>
template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_INTEGER(n_data);
  DATA_VECTOR(Y); 
  DATA_MATRIX(X); 
  DATA_MATRIX(B); // n_basis x n_data
  
  PARAMETER_VECTOR(ao); 
  PARAMETER_VECTOR(a); // length n_basis
  PARAMETER(ln_sig); 

  vector<Type> y_pred = X*ao + B*a; 

  Type obj = 0; 
  for(int i = 0; i < n_data; i++){
    obj -= dnorm(Y(i), y_pred(i), exp(ln_sig), true); 
  }

  return obj; 
}


