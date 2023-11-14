data {
    int<lower=0> Nrow;
    int<lower=0> Ncol;
    int<lower=0> Nplot; //Number of plots
    int<lower=0> Ncomp; // Number of non-missing values. 
    int<lower=0> Nmiss; // Number of missing values
    real dat_complete[Ncomp];   // Vector of non-missing values
    int ind_pres[Ncomp, 2];     // Matrix (row, col) of non-missing value indices
    int ind_miss[Nmiss, 2];     // Matrix (row, col) of missing value indices

    int<lower=0> Nrow_z;
    int<lower=0> Ncol_z;
    int<lower=0> Ncomp_z; // Number of non-missing values. 
    int<lower=0> Nmiss_z; // Number of missing values
    
    int idx_plot[Nrow];
    real dat_completez[Ncomp_z];   // Vector of non-missing values
    int ind_presz[Ncomp_z, 2];     // Matrix (row, col) of non-missing value indices
    int ind_missz[Nmiss_z, 2];     // Matrix (row, col) of missing value indices
    real tmaxAprMayJunscaled[Nrow, Ncol];
    real wateryrscaled[Nrow, Ncol];
    real SDI[Nrow];
    real MAP[Nplot];
    
}
parameters {
    // Multivariate normal distribution parameters
    real mu;
    real mu_plot;
    real<lower=0, upper =5> sigma_inc;
    real<lower=0, upper =5> sigma_add;
    real<lower=0, upper =5> sigma_dbh;
    // Vector containing "stochastic" nodes (for filling missing values
    real<lower=0> ymiss[Nmiss];
    real<lower=0> zmiss[Nmiss_z];
    real<lower=0> inc[Nrow, Ncol];
    real<lower=0> xinit[Nrow];
    real alpha_TREE[Nrow];
    real alpha_PLOT[Nplot];
    real<lower=1e-6> sigma_TREE;
    real<lower=1e-6> sigma_PLOT;
    
  
    real beta_YEAR[Ncol];
    real<lower=1e-6> sigma_YEAR;
 
    real betaTmax;
    real betaPrecip;
    real betaX;
    real betaSDI;
    real betaMAP;

    
  
   // real betaX_TREE[Nrow];

}

transformed parameters {
    real<lower=0> y[Nrow, Ncol];   // The "data" with interpolated missing values
    real<lower=0> z[Nrow, Ncol];   // The z "data" with interpolated missing values
    matrix<lower=0, upper = 125>[Nrow, Ncol] x;  // the estimated x values
   
    
     // Fill y with non-missing values 
    for(n in 1:Ncomp) {
        y[ind_pres[n,1], ind_pres[n,2]] = dat_complete[n];
    }
    // Fill the rest of y with missing value "parameters"
    for(n in 1:Nmiss){
        y[ind_miss[n, 1], ind_miss[n,2]] = ymiss[n];
    }
    
      // Fill x with non-missing values
     for(n in 1:Ncomp_z) {
         z[ind_presz[n,1], ind_presz[n,2]] = dat_completez[n];
     }
     // Fill the rest of y with missing value "parameters"
     for(n in 1:Nmiss_z){
         z[ind_missz[n, 1], ind_missz[n,2]] = zmiss[n];
     }
// estimate the true x value

 for(i in 1:Nrow){
   x[i,1] = xinit[i] + inc[i,1];
   
   for(t in 2:Ncol){ 
    x[i,t] = x[i,t-1] + inc[i,t];
    }
  }




  
}


model {
  
  for(i in 1:Nrow) {
    alpha_TREE[i] ~ normal(mu, sigma_TREE);
  }
  
  for(i in 1:Nplot) {
    alpha_PLOT[i] ~ normal(mu_plot, sigma_PLOT);
  }
 
// variance priors
 sigma_dbh ~ normal(1, 0.01); 
 sigma_inc ~ normal(0.035, 0.01);
 sigma_add ~ uniform(0, 5);
 sigma_PLOT ~ uniform(0, 5);
 
//fixed effect priors
 betaTmax ~ normal(0, 10);
 betaPrecip ~ normal(0, 10);
 betaX ~ normal(0, 10);
 betaSDI ~ normal(0, 10);
 mu_plot ~ normal(0, 10);
 


 
//data & process model for increment
  for(i in 1:Nrow){
    
      // initial condition for diameter
    xinit[i]  ~ uniform(0, 75);
      y[i,1] ~ normal(inc[i,1], sigma_inc)T[0,];
        inc[i,1] ~ normal(mu, sigma_add);
      z[i,1] ~ normal(x[i,1], 5); //or normal(x[i,1], sigma_dbh);

     for(t in 2:Ncol){
       //process model linear equation
       //Add in SDI effect here:
       inc[i,t] ~ lognormal(alpha_TREE[i] + alpha_PLOT[idx_plot[i]]+ betaX*x[i,t-1]+
       betaTmax*tmaxAprMayJunscaled[i,t]+ 
       betaPrecip*wateryrscaled[i,t] + betaSDI*SDI[i]+betaMAP*MAP[idx_plot[i]] , sigma_add);// + betaPrecip_MAT*wateryrscaled[i,t]*MAT[i] + betaSDI*SDI[i,t] + betaPrecip_Tmax*wateryrscaled[i,t]*tmaxAprMayJunscaled[i,t], sigma_add);
       //betaPrecip*wateryrscaled[idx_plot[i],t] //how you would index a time-varying climate matrix at the plot-level 
      
       //increment data model
         y[i,t] ~ normal(inc[i,t], sigma_inc)T[0,];
        //diameter data model
         z[i,t] ~ normal(x[i,t], sigma_dbh);
         
    }
  }
  
//print("xinit = ", xinit)
}
// generated quantities {
//   real log_lik[Nrow, Ncol];
//   for(i in 1:Nrow){
//     for(t in 1:Ncol) {
//       log_lik[i,t] = lognormal_lpdf(y[i,t] | inc[i,t] , sigma_add);
//     }
//   }
//}
