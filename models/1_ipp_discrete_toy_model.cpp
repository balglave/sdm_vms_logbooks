// Discrete version of a species distribution model accounting for preferential sampling
#include <TMB.hpp>

// Na values
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// square
template<class Type>
Type square(Type x){
  return pow(x,2.0);
}

// dlognorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=false){
  Type Return;
  if(give_log==false) Return = dnorm( log(x), meanlog, sdlog, false) / x;
  if(give_log==true) Return = dnorm( log(x), meanlog, sdlog, true) - log(x);
  return Return;
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  
  using namespace R_inla;
  using namespace Eigen;
  using namespace density;

  Type jnll=0;
  vector<Type> jnll_comp(4);
  jnll_comp.setZero();
  
  // Data
  DATA_FACTOR( Options_vec );
  DATA_INTEGER( n_cells );
  DATA_VECTOR( c_x );
  DATA_VECTOR( y_i );
  DATA_FACTOR( index_i );
  DATA_STRUCT( spde, spde_t );
  
  
  // Parameters
  PARAMETER( intercept_S );
  PARAMETER_VECTOR( deltainput_x );
  PARAMETER( intercept_l );
  PARAMETER( b );
  PARAMETER_VECTOR( etainput_x );
  PARAMETER_VECTOR( logtau );
  PARAMETER_VECTOR( logkappa );
  PARAMETER( logSD_obs );
  
  // Derived values
  int n_nodes = deltainput_x.size();
  int n_i = y_i.size();
  int n_gf = logkappa.size();
  Type SD_obs = exp(logSD_obs);
  vector<Type> MargSD(n_gf); MargSD.setZero();
  vector<Type> Range(n_gf); Range.setZero();
  
  for( int s=0; s<n_gf; s++){
    
    MargSD(s) = 1 / sqrt(4*M_PI) / exp(logtau(s)) / exp(logkappa(s));
    Range(s) = sqrt(8) / exp( logkappa(s) );
  
  }
  
  vector<Type> delta_x( n_nodes ); delta_x.setZero();
  delta_x = deltainput_x / exp(logtau(0)); // Transform random effects
  SparseMatrix<Type> Q_delta;
  Q_delta = Q_spde(spde, exp(logkappa(0)));
  if(Options_vec(2)==1) jnll_comp(0) = GMRF(Q_delta)(deltainput_x);
  
  vector<Type> eta_x( n_nodes ); eta_x.setZero();
  eta_x = etainput_x / exp(logtau(1)); // Transform random effects
  SparseMatrix<Type> Q_eta;
  Q_eta = Q_spde(spde, exp(logkappa(1)));
  if(Options_vec(2)==1) jnll_comp(1) = GMRF(Q_eta)(etainput_x);
  
  // Biomass field
  //--------------
  vector<Type> S_x(n_cells); S_x.setZero();
  for(int x=0; x<n_cells; x++){
    
    S_x(x) = exp(intercept_S + delta_x(x));
  
  }
  
  // Observation model
  //------------------
  for(int i=0; i<n_i; i++){
    if( !isNA(y_i(i)) & y_i(i) > 0 ){
      
      jnll_comp(2) -= dlognorm(y_i(i), log(S_x(index_i(i))) - square(SD_obs)/2, SD_obs, true);
    
    }
  }
  
  // Sampling process (Inhomogenous point process)
  //----------------------------------------------
  vector<Type> log_lambda_x(n_cells);log_lambda_x.setZero();
  vector<Type> lambda_x(n_cells);
  lambda_x = exp(log_lambda_x);
  vector<Type> loglambda_i(n_i);loglambda_i.setZero();
  
  // Raw form
  if(Options_vec(0)==0){
    for(int i=0; i<n_i; i++){

      // likelihood of data points
      loglambda_i(i) = intercept_l + b * log(S_x(index_i(i))) + eta_x(index_i(i));
      if(Options_vec(3)==1) jnll_comp(3) -= loglambda_i(i);

    }
    
    // integrating constant of the ipp
    for(int x=0; x<n_cells; x++){
      
      lambda_x(x) = exp( intercept_l + b * log(S_x(x)) + eta_x(x) );
      if(Options_vec(3)==1) jnll_comp(3) -= (Type(1)-lambda_x(x));
      
    }
    
  }
  
  // Conditionnal form
  if(Options_vec(0)==1){
    for(int x=0; x<n_cells; x++){
      log_lambda_x(x) =  intercept_l + b * log(S_x(x)) + eta_x(x) ;
      if(!isNA(c_x(x)) & Options_vec(3)==1) jnll_comp(3) -= log_lambda_x(x) - log(lambda_x.sum());
    }
  }
  
  // Report
  REPORT(S_x);
  REPORT(intercept_S);
  REPORT(delta_x);
  REPORT(lambda_x);
  REPORT(intercept_l);
  REPORT(b);
  REPORT(eta_x);
  REPORT(Range);
  REPORT(MargSD);
  
  if(Options_vec(1)==1){
    ADREPORT(S_x);
  }
  
  // sum over likelihood components
  int n_jnll = jnll_comp.size();
  for(int j=0; j<n_jnll; j++){
    jnll += jnll_comp(j);
  }
  
  return jnll;

}
