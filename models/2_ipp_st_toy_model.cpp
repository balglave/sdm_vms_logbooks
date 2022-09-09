// Spatio-temporal species distribution model accounting for preferential sampling
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
  
  // Negative log-likelihood
  Type jnll=0;
  vector<Type> jnll_comp(4);
  jnll_comp.setZero();
  
  // Model inputs
  //-------------
  // Data
  DATA_FACTOR( Options_vec );
  DATA_INTEGER( n_cells );
  DATA_INTEGER( n_t );
  DATA_INTEGER( n_i );
  DATA_INTEGER( n_gf );
  DATA_INTEGER( n_nodes );
  DATA_VECTOR( y_i );
  DATA_FACTOR( t_i );
  DATA_IMATRIX( Aix_ij );
  DATA_VECTOR( Aix_w );
  DATA_IMATRIX( Aix_ij_pred );
  DATA_VECTOR( Aix_w_pred );
  DATA_STRUCT( spde, spde_t );
  
  // Parameters
  PARAMETER_VECTOR( intercept_S );
  PARAMETER_MATRIX( deltainput_x );
  PARAMETER_VECTOR( intercept_l );
  PARAMETER_VECTOR( b );
  PARAMETER_MATRIX( etainput_x );
  PARAMETER_VECTOR( logtau );
  PARAMETER_VECTOR( logkappa );
  PARAMETER( logSD_obs );
  PARAMETER( rho_delta );
  
  // Derived values
  int i,x;
  Type SD_obs = exp(logSD_obs);
  vector<Type> MargSD(n_gf);
  vector<Type> Range(n_gf);
  
  for( int s=0; s<n_gf; s++){
    
    MargSD(s) = 1 / sqrt(4*M_PI) / exp(logtau(s)) / exp(logkappa(s));
    Range(s) = sqrt(8) / exp( logkappa(s) );
  
  }
  
  
  // Random effects
  //---------------
  matrix<Type> delta_x( n_nodes, n_t );
  SparseMatrix<Type> Q_delta;
  
  matrix<Type> eta_x( n_nodes, n_t );
  SparseMatrix<Type> Q_eta;
  
  for( int t=0; t<n_t; t++){
    
    // Transform random effects
    for( int x=0; x<n_nodes; x++){
      
      delta_x(x,t) = deltainput_x(x,t) / exp(logtau(0));
      eta_x(x,t) = etainput_x(x,t) / exp(logtau(1));
      
    }
    
    Q_delta = Q_spde(spde, exp(logkappa(0)));
    Q_eta = Q_spde(spde, exp(logkappa(1)));
    
    if(t == 0) jnll_comp(0) += GMRF(Q_delta)(deltainput_x.col(t));
    if(t > 0) jnll_comp(0) += GMRF(Q_delta)(deltainput_x.col(t) - rho_delta*deltainput_x.col(t-1));
    
    jnll_comp(1) +=  GMRF(Q_eta)(etainput_x.col(t));
  
  }
  
  // Link to data
  //-------------
  vector<Type> delta_i(n_i); delta_i.setZero();
  vector<Type> eta_i(n_i); eta_i.setZero();
  for( int Arow=0; Arow<Aix_ij.rows(); Arow++ ){
    
    i = Aix_ij(Arow,0);
    x = Aix_ij(Arow,1);
    delta_i(i) += Aix_w(Arow) * delta_x(x,t_i(i));
    eta_i(i) += Aix_w(Arow) * eta_x(x,t_i(i));
    
  }
  
  // Likelihood of the observations
  //-------------------------------
  vector<Type> S_i(n_i); S_i.setZero();
  vector<Type> log_lambda_i(n_i); log_lambda_i.setZero();
  for(int i=0; i<n_i; i++){
    
    S_i(i) = exp(intercept_S(t_i(i)) + delta_i(i));
    
    // catch
    if( !isNA(y_i(i)) & y_i(i) > 0 ){
      
      jnll_comp(2) -= dlognorm(y_i(i), log(S_i(i)) - square(SD_obs)/2, SD_obs, true);
      
    }
    
    // sampling
    log_lambda_i(i) = intercept_l(t_i(i)) + b(t_i(i)) * log(S_i(i)) + eta_i(i);
    jnll_comp(2) -=  log_lambda_i(i); // log intensity of each data point ( see expression of an ipp in Diggle (2013))
    // n.b the integrating constant (integration over the full area) of the ipp
    // is computed at the scale of the prediction grid (see below when 
    // projecting fishing intensity on the prediction grid)
  }
  
  // Projecting biomass field and sampling intensity 
  // field over the projection grid
  //------------------------------------------------
  matrix<Type> delta_p( n_cells, n_t ); delta_p.setZero();
  matrix<Type> eta_p( n_cells, n_t ); eta_p.setZero();
  matrix<Type> S_p(n_cells,n_t); S_p.setZero();
  matrix<Type> lambda_p(n_cells,n_t); lambda_p.setZero();
  for(int t=0; t<n_t; t++){
    
    // Projection of random effect on the prediction grid
    for( int Arow=0; Arow<Aix_ij_pred.rows(); Arow++ ){
      
      int p = Aix_ij_pred(Arow,0);
      int x = Aix_ij_pred(Arow,1);
      delta_p(p,t) += Aix_w_pred(Arow) * delta_x(x,t);
      eta_p(p,t) += Aix_w_pred(Arow) * eta_x(x,t);
      
    }
    
    // Projection of biomass field
    for(int p=0; p<n_cells; p++){

      S_p(p,t) = exp(intercept_S(t) + delta_p(p,t) );
      lambda_p(p,t) = exp(intercept_l(t) + b(t) * log(S_p(p,t)) + eta_p(p,t) );
      jnll_comp(2) -= Type(1)-lambda_p(p,t); // integrating constant
      
    }
  }
  
  // Report
  REPORT(S_p);
  REPORT(intercept_S);
  REPORT(delta_x);
  REPORT(delta_p);
  REPORT(lambda_p);
  REPORT(intercept_l);
  REPORT(b);
  REPORT(eta_x);
  REPORT(eta_p);
  REPORT(Range);
  REPORT(MargSD);
  REPORT(rho_delta);
  REPORT(SD_obs);
  
  if(Options_vec(0)==1){
    ADREPORT(S_p);
  }
  
  
  // sum over likelihood components
  int n_jnll = jnll_comp.size();
  for(int j=0; j<n_jnll; j++){
    jnll += jnll_comp(j);
  }
  
  return jnll;

}
