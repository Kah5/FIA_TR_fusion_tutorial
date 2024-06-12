data {
  // ...
  int<lower=0> Nrow; //Ntrees
  int<lower=0> Ncol; //Nyears
  int<lower=0> Ndia; //total number of diameter measurements
  int<lower=0> Ninc; //total number of increment measurements
  int<lower=0> treedia[Ndia]; //index of z diameter measurements denoting the tree number
  int<lower=0> yeardia[Ndia]; //index of z diameter measurements denoting the year number
  real<lower = 0> z[Ndia];//diameter measurement vector

  real<lower = 0> y[Ninc]; // all the non-missing tree ring increments
  int<lower = 0> treeinc[Ninc];//index of y increment measurements denoting the tree number
  int<lower = 0> yearinc[Ninc];//index of y increment measurements denoting the year number
  matrix [Nrow, Ncol] tmaxAprMayJunscaled;
  matrix [Nrow, Ncol] wateryrscaled;
  // ...
}
parameters {
    // variance parameters
  
    real<lower=0, upper =10> sigma_inc;
    real<lower=0, upper =10> sigma_add;
    real<lower=0, upper =10> sigma_dbh;
    
    //tree level random effect parameters
    real<lower=-10, upper =10> mutree;
    real alpha_TREE[Nrow];
    real<lower=1e-6> sigma_TREE;
    
   
   //setting inc and xinit as parameters
    matrix<lower=0, upper =2> [Nrow, Ncol] inc;
    vector<lower=0>[Nrow] xinit;
  
    
    real betaX;
    real betaTmax;
    real betaPrecip;


}

transformed parameters{
   matrix<lower=1e-6>[Nrow, Ncol] x; //true diameter x is a transformed parameter
   matrix[Nrow, Ncol+1] xscaled;

  for (i in 1:Nrow){
    //x[i,1] = xinit[i] + incinit[i];
    x[i,1] = xinit[i] + inc[i, 1];

    for (t in 2:Ncol) {
      x[i, t] = x[i, t-1] + inc[i, t] ;
    }
}

//rescaling diameter helps with estimation and speed of betaX
//also putting the initial diameter and the estimated diameters in the same matrix helps alot
//I think this speeds things up by only having a single line defining the liklihood function
xscaled[1:Nrow, 1] = (xinit-31.63)/10.61;
xscaled[,2:(Ncol+1)] = (x - 31.63)/10.61; //scale diameter by the mean and sd 
}
model {
  
  //priors on tree-level random effect
  mutree ~ normal(0, 5);
  sigma_TREE ~ cauchy(0, 5);
 
  alpha_TREE[1:Nrow] ~ normal(mutree, sigma_TREE);
 
  
  //variance priors
  sigma_dbh ~ normal(1, 0.01); //normal(1,0.01) works well on base model
  sigma_inc ~ normal(0.035, 0.01);
  sigma_add ~ normal(0, 5);

  betaX ~normal(0,10);
  betaTmax ~ normal(0, 10);
  betaPrecip ~ normal(0, 10);
 
 
 // initial condition prior for diameter, for each tree
      
  xinit  ~ uniform(0, 75);
  
       
  
 // tree ring diameter increment process model
 for(i in 1:Nrow){//loop over the number of trees
     
        inc[i,1:Ncol] ~ lognormal( alpha_TREE[i] + betaX*xscaled[i, 1:(Ncol)] + betaTmax*tmaxAprMayJunscaled[i,1:Ncol] + betaPrecip*wateryrscaled[i,1:Ncol] , sigma_add);
      
  
 }
  //tree ring increment data model
  for(j in 1:Ninc){ //loop over the number of increment measurements
     y[j] ~ normal(inc[treeinc[j], yearinc[j]], sigma_inc)T[0,];
     //increment data y corresponds to the true increment matrix through the tree and year increment indices
   }
 
 //diameter data model --moved from previous loop
   //  z[i,t] ~ normal(x[i,t], sigma_dbh);
   for(d in 1:Ndia){
     z[d] ~ normal(x[treedia[d], yeardia[d]], sigma_dbh);
     //diameter data z corresponds to the true diameter matrix through the tree and year diameter indices
   }
}
