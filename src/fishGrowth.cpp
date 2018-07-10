#include <TMB.hpp> 
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{ 
  DATA_INTEGER(ndiv); //number of divisions
	DATA_INTEGER(nobs); //number of observations
	DATA_VECTOR(Lobs);     //length observations
	DATA_VECTOR(age); 
	DATA_VECTOR(Wobs); //weight observations
	DATA_IVECTOR(divs); //index for division
	DATA_INTEGER(nc); //number of cohorts
	DATA_IVECTOR(ic); //index for cohorts
	//DATA_VECTOR(log_L0);//length at hatch
//Covariates
	DATA_VECTOR(sx); //index for sex - 0 female, 1 male
	DATA_VECTOR(mat); //index for maturity - 0 immature, 1 mature
	
//For year-effects	
	//DATA_INTEGER(ny); //number of year effects
	//DATA_IVECTOR(iy); //index for year effects
	
	
	
	
	
//Fixed effects
   PARAMETER(log_Linf); 
   PARAMETER(log_k); 
   PARAMETER(t0);    
   PARAMETER(log_a);  
   PARAMETER(log_b);  
   PARAMETER(log_std_me1); //observation error length
   PARAMETER(log_std_me2); //observation error weight
   //std deviation between divisions
   PARAMETER(log_std_log_Linf); 
   PARAMETER(log_std_log_k);
   PARAMETER(log_std_log_a);
   PARAMETER(log_std_log_b);
   //std deviation between cohorts
   PARAMETER_VECTOR(log_std_log_Linfc);
   PARAMETER_VECTOR(log_std_log_kc);
   PARAMETER_VECTOR(log_std_log_ac);
   PARAMETER_VECTOR(log_std_log_bc);
   
//Covariate effects
   PARAMETER(log_maty);
   PARAMETER(log_sex);
	
//Random Effects;
//division effects
   PARAMETER_VECTOR(rlog_Linf);
   PARAMETER_VECTOR(rlog_k);
   PARAMETER_VECTOR(rlog_a);
   PARAMETER_VECTOR(rlog_b);
//cohort effects   
   PARAMETER_MATRIX(rlog_Linfc);
   PARAMETER_MATRIX(rlog_kc);
   PARAMETER_MATRIX(rlog_ac);
   PARAMETER_MATRIX(rlog_bc);
   
//For year-effects	
	//PARAMETER(log_std_log_Linfy);
	//PARAMETER_VECTOR(rlog_Linfy);
	//Type std_log_Linfy = exp(log_std_log_Linfy);

//Other declarations for conversion from log space etc	
//division effect standard deviation
   Type std_log_Linf = exp(log_std_log_Linf);
   Type std_log_k = exp(log_std_log_k);
   Type std_log_a = exp(log_std_log_a);
   Type std_log_b = exp(log_std_log_b);

//cohort effect standard deviation
   vector<Type> std_log_Linfc = exp(log_std_log_Linfc);
   vector<Type> std_log_kc = exp(log_std_log_kc);
   vector<Type> std_log_ac = exp(log_std_log_ac);
   vector<Type> std_log_bc = exp(log_std_log_bc);

//observation error standard deviation
   Type std_me1 = exp(log_std_me1); 
   Type std_me2 = exp(log_std_me2); 
  
 //k  
 //  Type k = exp(log_k); 

  
  //Declarations
  //divisions 
   vector<Type> log_Linf_div = log_Linf + rlog_Linf;
   vector<Type> log_k_div = log_k + rlog_k;
   vector<Type> log_a_div = log_a + rlog_a;
   vector<Type> log_b_div = log_b + rlog_b;
  //divisions and cohorts 
   array<Type> log_Linf_div_c(ndiv,nc);
   array<Type> log_k_div_c(ndiv,nc);
   array<Type> log_a_div_c(ndiv,nc);
   array<Type> log_b_div_c(ndiv,nc);
   array<Type> log_Winf_div_c(ndiv,nc);
   array<Type> Winf_div_c(ndiv,nc);
   array<Type> Linf_div_c(ndiv,nc);
   array<Type> k_div_c(ndiv,nc);

   //filling in arrays
   for(int i = 0;i < ndiv;++i){
     for(int j = 0;j < nc;++j)
     {log_Linf_div_c(i,j) = log_Linf_div(i) + rlog_Linfc(i,j);
	  log_k_div_c(i,j) = log_k_div(i) + rlog_kc(i,j);
	  log_a_div_c(i,j) = log_a_div(i) + rlog_ac(i,j);
	  log_b_div_c(i,j) = log_b_div(i) + rlog_bc(i,j);
	  log_Winf_div_c(i,j) =log_a_div_c(i,j)+exp(log_b_div_c(i,j))*log_Linf_div_c(i,j);
      Linf_div_c(i,j) = exp(log_Linf_div_c(i,j));
      k_div_c(i,j) =exp(log_k_div_c(i,j));
	  Winf_div_c(i,j) =exp(log_Winf_div_c(i,j));	 }  }
   
 
   
   vector<Type> log_Lobs = log(Lobs);   
   vector<Type> log_Wobs = log(Wobs); 
   
   vector<Type> Lhat(nobs);   
   vector<Type> log_Lhat(nobs);
   vector<Type> resid1(nobs);  
   vector<Type> std_resid1(nobs);

   vector<Type> What(nobs);   
   vector<Type> log_What(nobs);
   vector<Type> resid2(nobs);  
   vector<Type> std_resid2(nobs);
   vector<Type> log_linf_all(nobs);
   vector<Type> log_k_all(nobs);
   vector<Type> k_all(nobs);
   vector<Type> log_a_all(nobs);
   vector<Type> b_all(nobs);
   
//Applying Covariate effects on k or b
   vector<Type> sex_all;
   sex_all=exp(log_sex)*sx;

   vector<Type> mat_all;
   mat_all=exp(log_maty)*mat;
   

  
   Type zero = 0.0;
   //Type one = 1.0; //for the Winf model
   
//matching with division and cohort in data
   Type nll = zero;
   for(int j = 0;j < nobs;++j)
   {
   log_linf_all = log_Linf_div_c(divs(j),ic(j));
   log_k_all(j) = log_k_div_c(divs(j),ic(j))+sex_all(j)+mat_all(j);
   k_all(j)=exp(log_k_all(j));
   log_a_all(j) = log_a_div_c(divs(j),ic(j));
   b_all(j) = exp(log_b_div_c(divs(j),ic(j)))+mat_all(j);}
   
 //vbgf model with t0  
   log_Lhat = log_linf_all + log(1-exp(-k_all*(age-t0))); 
 //model providing without t0 - length of hatchlings
 //log_Lhat = log_L0-(k_all*age)+log_linf_all + log(1-exp(-k_all*age)); 

   Lhat = exp(log_Lhat);
   resid1 = log_Lobs - log_Lhat;
   std_resid1 = resid1/std_me1;    
   
   //log_What = log_a_all + b_all*log_Lhat; //uses predicted lengths
   log_What = log_a_all + b_all*log_Lobs; //uses observed lengths
   What = exp(log_What);
   resid2 = log_Wobs - log_What;
   std_resid2 = resid2/std_me2;  

//model with weight vonb - when weight at age data is available
//needs log_Winf and po as parameters 
//is it possible to use existing declarations of log_Linf and t0?
   //log_What = log_Winf_all + b_all*log(one - (one - po)*exp(-k_all*age));    
   
// nll for data;
   nll -= dnorm(resid1, zero, std_me1, true).sum();
   nll -= dnorm(resid2, zero, std_me2, true).sum();

// nll for division effects   
   nll -= dnorm(rlog_Linf, zero, std_log_Linf, true).sum(); 
   nll -= dnorm(rlog_k, zero, std_log_k, true).sum(); 
   nll -= dnorm(rlog_a, zero, std_log_a, true).sum(); 
   nll -= dnorm(rlog_b, zero, std_log_b, true).sum(); 
   
//nll for Linf cohort effects   
   for (int i = 0; i < ndiv; ++i){
   nll -= dnorm(rlog_Linfc(i,0), zero, std_log_Linfc(i), true); 
   for (int j = 1; j < nc; ++j) {  
     nll -= dnorm(rlog_Linfc(i,j), rlog_Linfc(i,j-1), std_log_Linfc(i), true);
   } }

 //nll for k cohort effects   
   for (int i = 0; i < ndiv; ++i){
   nll -= dnorm(rlog_kc(i,0), zero, std_log_kc(i), true); 
   for (int j = 1; j < nc; ++j) {  
     nll -= dnorm(rlog_kc(i,j), rlog_kc(i,j-1), std_log_kc(i), true);
   } }
  
   
   //nll for log_a cohort effects   
   for (int i = 0; i < ndiv; ++i){
   nll -= dnorm(rlog_ac(0), zero, std_log_ac(i), true); 
   for (int j = 1; j < nc; ++j) {  
     nll -= dnorm(rlog_ac(j), rlog_ac(j-1), std_log_ac(i), true);
   } }
 
//nll for b cohort effects   
   for (int i = 0; i < ndiv; ++i){
   nll -= dnorm(rlog_bc(0), zero, std_log_bc(i), true); 
   for (int j = 1; j < nc; ++j) {  
     nll -= dnorm(rlog_bc(j), rlog_bc(j-1), std_log_bc(i), true);
   } }
   



   REPORT(Lhat);
   REPORT(log_Lhat);
   REPORT(resid1);
   REPORT(std_resid1);
   REPORT(resid2);
   REPORT(std_resid2);
   REPORT(k_div_c);  
   REPORT(Winf_div_c);
   REPORT(Linf_div_c);
   REPORT(log_a_div_c);
   REPORT(log_b_div_c);
   
   ADREPORT(log_Linf_div_c);
   ADREPORT(log_k_div_c);
   ADREPORT(log_a_div_c);
   ADREPORT(log_b_div_c);
   ADREPORT(log_Winf_div_c);


   return nll;
}















