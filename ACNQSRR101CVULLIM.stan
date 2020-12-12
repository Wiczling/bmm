functions{
  real hplcmodel(real fi, real logkw, real logka, real logS2){
    
    real logk;												// retention factor
    real S1;
    
    S1 = (logkw - logka)*(1+10^logS2);
    logk = logkw - S1 * fi / (1 + 10^logS2 * fi);
    
    return logk;
  }
}

data{
  int nAnalytes;			// number of analytes
  int nObs;					// number of observations
  int analyte[nObs];		// analytes indexes
  int start[nAnalytes];		// first apperance of analyte in "analyte" vector
  vector[nObs] logkObs;	    // observed retention factors
  vector[nObs] fi;			// organic modifier content in the mobile phase
  real Mmolx[nObs];         // (moleculuar mass-300)/100
  int<lower=0> K;           //  number of predictors (functional groups)
  matrix[nAnalytes, K] nrfungroups;   // predictor matrix (functional groups)   
  int<lower = 0, upper = 1> run_estimation; // 0 for prior predictive, 1 for estimation
  int nEst;			       // number of observations
  int cvidx[nEst];		   // analytes indexes
  vector[nEst] logkObsEst; // observed retention factors
  int nEstmax;			       // number of observations
  int cvidxmax[nEstmax];		   // analytes indexes
  int nEstmin;			       // number of observations
  int cvidxmin[nEstmin];		   // analytes indexes
  real fix;
}

parameters{

 real logkwHat;	// mean value of logkw 
 real logkaHat; // mean value of logka 
 real logS2Hat;	// mean curvature coefficient 
 real<lower = 0> sigma;		// standard deviation for residuals
 vector<lower = 0>[3] omega;// diagonal elements of variance-covariance matrix for inter-analyte variability 
 corr_matrix[3] rho;	    // correlation matrix		
 real<lower = 1> nu;	    // normality constant for inter-analyte variability 
 real<lower = 1> nuobs;     // normality constant for residual variability 
 real<lower = 1> nupi;     // normality constant for residual variability 
 real beta[3];			    // regression coefficients for Mmolx
 vector[3] param[nAnalytes]; // individual values of chromatographic parameters
 vector<lower = 0>[K] pilogkw;  // regression coefficient for logkw
 vector[K] pidlogk ;  //... logka logkw difference
 vector[K] pilogS2;             // ... logS2
 real<lower = 0> spilogkw;      // group-level std for logkw
 real<lower = 0> spidlogk;      //... logka
 real<lower = 0> spilogS2;      //... logS2
 real<lower = 0> mpilogkw;      // group-level mean for logkw
 real mpidlogk;      //... logka
}

transformed parameters{
  vector[3] miu[nAnalytes];	 
  real logka[nAnalytes];
  real logkw[nAnalytes];
  real logS2[nAnalytes];
  vector[K] pilogka;
  cov_matrix[3] Omega;			 // variance-covariance matrix
  vector[nObs] logkHat;		
  vector[nEst] logkHatEst;
  vector[nEstmax] logkHatEstmax;	
  vector[nEstmin] logkHatEstmin;	


  Omega = quad_form_diag(rho, omega);	// diag_matrix(omega) * rho * diag_matrix(omega)

  pilogka = pilogkw - pidlogk;

  for(j in 1:nAnalytes){
    miu[j,1]  = logkwHat + beta[1] * Mmolx[start[j]] - nrfungroups[j,1:K] * pilogkw;
    miu[j,2]  = logkaHat + beta[2] * Mmolx[start[j]] - nrfungroups[j,1:K] * pilogka; 
    miu[j,3]  = logS2Hat + beta[3] * Mmolx[start[j]] + nrfungroups[j,1:K] * pilogS2;
  }

	for(j in 1:nAnalytes){
		logkw[j] = param[j, 1];
		logka[j] = param[j, 2];
		logS2[j] = param[j, 3];
	}
  
  for(i in 1:nObs){
    logkHat[i] = hplcmodel(fi[i], logkw[analyte[i]], logka[analyte[i]], logS2[analyte[i]]);
 }

 logkHatEst = logkHat[cvidx];

 for(i in 1:nEstmin){
 logkHatEstmin[i] =  hplcmodel(fix, logkw[analyte[cvidxmin[i]]], logka[analyte[cvidxmin[i]]], logS2[analyte[cvidxmin[i]]]);
 }

 for(i in 1:nEstmax){
 logkHatEstmax[i] =  hplcmodel(fix, logkw[analyte[cvidxmax[i]]], logka[analyte[cvidxmax[i]]], logS2[analyte[cvidxmax[i]]]);
 }

}

model{
logkwHat ~ normal(6.6, 1.5);  //3.6+2*1.5
  logkaHat ~ normal(1.3, 1.5); //-1.7+2*1.5
  logS2Hat ~ normal(log10(2), 0.2);

  beta[1]  ~ normal(1.4,1.5);
  beta[2]  ~ normal(0.2,1.5);
  beta[3]  ~ normal(0,0.2);

  omega[1] ~ normal(0,1.50);
  omega[2] ~ normal(0,1.50);
  omega[3] ~ normal(0,0.2);

  rho   ~ lkj_corr(1);
  sigma  ~ normal(0,0.067);

  mpilogkw ~ normal(0,1.5);
  mpidlogk ~ normal(0,1.5);
  
  spilogkw ~ normal(0,1.5);
  spidlogk ~ normal(0,1.5);
  spilogS2 ~ normal(0,0.2);

  pilogkw ~ lognormal(log(mpilogkw),spilogkw);
  pidlogk ~ student_t(nupi,mpidlogk,spidlogk);
  pilogS2 ~ normal(0,spilogS2);

  nu    ~ gamma(2,0.1);
  nuobs ~ gamma(2,0.1);
  nupi  ~ gamma(2,0.1);

  for(i in  1:nAnalytes){
   param[i] ~ multi_student_t(nu,miu[i],Omega);
  }
  
  if(run_estimation==1){
   logkObsEst ~ student_t(nuobs,logkHatEst, sigma);// likelihood
   target += student_t_lccdf(2|nuobs,logkHatEstmax, sigma);  ## censored data likelihood
   target += student_t_lcdf(-1.5|nuobs,logkHatEstmin, sigma);  ## censored data likelihood
  }
}

generated quantities{
  real logkCond[nObs];
  real logkPred[nObs];
  //real log_lik[nObs];
  vector[3] paramPred[nAnalytes]; 
  
  for(j in 1:nAnalytes){
   paramPred[j] = multi_student_t_rng(nu,miu[j],Omega);
  }
  
  for(i in 1:nObs){
   real logkHatPred;	// predicted logk	
   logkHatPred = hplcmodel(fi[i], paramPred[analyte[i],1], paramPred[analyte[i],2], paramPred[analyte[i],3]);
   logkCond[i] = student_t_rng(nuobs, logkHat[i], sigma);
   logkPred[i] = student_t_rng(nuobs, logkHatPred, sigma);
   //log_lik[i]= student_t_lpdf(logkObs[i] | nuobs, logkHat[i], sigma);
  }
}
