#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector normalize_vector(NumericVector p){
  int n = p.size();
  double sump = sum(p);
  double inv_sum_p = 1/sump;
  NumericVector out = p;
  
  for(int i = 0; i < n; i++){
    out[i] *= inv_sum_p; 
  }
  
  return out;
}


// [[Rcpp::export]]
NumericMatrix compute_joint_prob(NumericVector probf,
                                 NumericVector probm, 
                                 NumericVector corr){
  
  int k = probf.size();
  NumericMatrix prob_matrix(4.0, k);
  
  for(int t = 0; t < k; t++){
    double sigF = pow(probm[t] * (1-probm[t]), 0.5);
    double sigM = pow(probf[t] * (1-probf[t]), 0.5);
    double corr_effect = corr[t] * sigF * sigM;
    prob_matrix(0, t) = (1-probf[t]) * (1-probm[t]) + corr_effect;
    prob_matrix(1, t) = probf[t] * (1-probm[t])     - corr_effect;
    prob_matrix(2, t) = (1-probf[t]) * probm[t]     - corr_effect;
    prob_matrix(3, t) = probf[t] * probm[t]         + corr_effect;
  }
  
  return prob_matrix;
}

// [[Rcpp::export]]
NumericMatrix compute_cond_prob(NumericMatrix JointProb,
                                NumericVector probf,
                                NumericVector probm){

  int k = probf.size();
  NumericMatrix prob_matrix(2.0, k);

  for(int t = 0; t < k; t++){
    prob_matrix(0, t) = JointProb(1,t)/(1-probm[t]);
    prob_matrix(1, t) = JointProb(3,t)/probm[t];
  }

  return prob_matrix;
}


double minIntCpp(int x,
                 int y){
  if(x > y){
    return(y);
  } else {
    return(x);
  }
  
}

double maxdblCpp(double x,
                 double y){
  if(x > y){
    return(x);
  } else {
    return(y);
  }
  
}


double maxIntCpp(int x,
                 int y){
  if(x > y){
    return(x);
  } else {
    return(y);
  }
  
}


// [[Rcpp::export]]
double minIntCppV(IntegerVector x){
  
  int n = x.size();
  double out;
  out = x[0];
  
  for(int i = 1; i < n; i++){
      out = minIntCpp(x[i],out); 
  }
  
  return out;
}

// [[Rcpp::export]]
IntegerVector initialize_first_capture(int n,
                         int k,
                         NumericVector tau){
  
  IntegerVector out = sample(k, n, true, tau); 
  
  return out;
}

NumericMatrix initialize_recruit(int n,
                                 int k,
                                 IntegerVector first_capture){
  
  NumericMatrix out(n,k);

  for(int i = 0; i < n; i++){
    int t0 = first_capture[i];
    for(int t = t0; t < k; t++){
      out(i,t) = 1;
    }
  }
  
  return out;  
}

NumericMatrix initialize_mating(int n,
                                int k,
                                IntegerVector first_capture,
                                NumericVector delta){
  
  NumericMatrix out(n,k);
  for(int i = 0; i < n; i++){
    int t0 = first_capture[i];
    out(i,t0) = rbinom(1,1,delta[t0])[0];
  }
  
  return out;
}
  
// [[Rcpp::export]]
NumericMatrix initialize_survival(int n,
                                  int k,
                                  IntegerVector first_capture){
  
  NumericMatrix out(n,k);
  for(int i = 0; i < n; i++){
    int t0 = first_capture[i];
    for(int t=0; t < (t0+1); t++){
      out(i,t) = 1;
    }
  }
  return out;
}

// [[Rcpp::export]]
IntegerMatrix initialize_pairs(int nf,
                               int nm,
                               int k,
                               int t0,
                               NumericMatrix mating_f,
                               NumericMatrix mating_m){
  
  IntegerMatrix out(nf,k);
  
  for(int i = 0; i < nf; i++){
    for(int t = 0; t < k; t++){
      out(i,t) = nm;
    }
  }
  
  NumericMatrix prob_mate(nf, nm);
  double females_mating = 0;
  double males_mating = 0; 
  
  // Probability that a pair can form 
  for(int i = 0; i < nf; i++){
    for(int j = 0; j < nm; j++){
      prob_mate(i,j) = mating_f(i,t0) * mating_m(j,t0);
    }
  }
  
  //Number of mating females
  for(int i = 0; i < nf; i++){
    if(mating_f(i,t0) == 1){
      females_mating += 1;
    }
  }

  //Number of mating males
  for(int j = 0; j < nm; j++){
    if(mating_m(j,t0) == 1){
      males_mating += 1;
    }
  }

  // Add a min of zero
  males_mating = maxdblCpp(males_mating, 1.0);

  // Operating Sex Ratio (which side randomly samples from bigger side)
  double osr = females_mating/males_mating;
  //Sample from larger pop in 2D without replacement
  if(osr <= 1){
    for(int i = 0; i < nf; i++){
      NumericVector prob = prob_mate(i, _ );
      double sump = sum(prob);
      
      // If no probability then just assign dummy state
      if(sump > 0){
        // prob = normalize_vector(prob);
        int j = sample(nm, 1, false, prob)[0]-1;
        out(i,t0) = j;
        prob_mate(i, _ ) = prob_mate(i, _ ) *  0;
        
        if(j != nm){
          prob_mate( _ ,j) = prob_mate( _ ,j) * 0;
        }
      }
      
    }
  } else {
    for(int j = 0; j < nm; j++){
      NumericVector prob = prob_mate( _ ,j);
      double sump = sum(prob);
      
      if(sump != 0){
        int i = sample(nf, 1, false, prob)[0]-1;
        out(i,t0) = j;
        
        if(i != nf){
          prob_mate(i, _ ) = prob_mate(i, _ ) *  0;
        }
        
        prob_mate( _ ,j) = prob_mate( _ ,j) * 0;
      }
    }
  }

  return out;
}



// [[Rcpp::export]]
IntegerVector compute_pairs_t(NumericVector repartner,
                              IntegerVector pairs_previous,
                              NumericVector mating_f,
                              NumericVector mating_m,
                              int nf,
                              int nm){
  
  IntegerVector out(nf, nm);
  NumericMatrix prob_mate(nf, nm);
  double females_mating = 0;
  double males_mating = 0; 
  
  //Number of mating females
  for(int i = 0; i < nf; i++){
    if(mating_f[i] == 1){
      females_mating += 1;
    }
  }
  
  //Number of mating males
  for(int j = 0; j < nm; j++){
    if(mating_m[j] == 1){
      males_mating += 1;
    }
  }
  
  // Add a min of zero
  males_mating = maxdblCpp(males_mating, 1.0);
  
  // Operating Sex Ratio (which side randomly samples from bigger side)
  double osr = females_mating/males_mating;
  
  // Probability that a pair can form 
  for(int i = 0; i < nf; i++){
    if(repartner[i] == 1.0){
      int j = pairs_previous[i];
      prob_mate(_, j) = prob_mate(_, j) * 0;
      prob_mate(i, _) = prob_mate(i, _) * 0;
      prob_mate(i,j) = 1;
      mating_f[i] = 0;
      mating_m[j] = 0;
    } else {
      for(int j = 0; j < nm; j++){
        prob_mate(i,j) = mating_f[i] * mating_m[j];
      }
    }
  }

  //Sample from larger pop in 2D without replacement
  if(osr <= 1){
    for(int i = 0; i < nf; i++){
      NumericVector prob = prob_mate(i, _ );
      double sump = sum(prob);

      // If no probability then just assign dummy state
      if(sump > 0){
        // prob = normalize_vector(prob);
        int j = sample(nm, 1, false, prob)[0]-1;
        out[i] = j;
        prob_mate(i, _ ) = prob_mate(i, _ ) *  0;

        if(j != nm){
          prob_mate( _ ,j) = prob_mate( _ ,j) * 0;
        }
      }

    }
  } else {
    for(int j = 0; j < nm; j++){
      NumericVector prob = prob_mate( _ ,j);
      double sump = sum(prob);

      if(sump > 0){
        int i = sample(nf, 1, false, prob)[0]-1;
        out[i] = j;

        if(i != nf){
          prob_mate(i, _ ) = prob_mate(i, _ ) *  0;
        }

        prob_mate( _ ,j) = prob_mate( _ ,j) * 0;
      }
    }
  }
  
  return out;
}


// [[Rcpp::export]]
NumericVector compute_mating_t(NumericVector surv,
                               IntegerVector first,
                               int t,
                               double delta){
  
  int n = first.size();
  NumericVector out(n);
  
  for(int i = 0; i < n; i++){
    if(first[i] < t){
      if(surv[i] == 1){
        out[i] = rbinom(1,1,delta)[0];
      }
    }
  }
  
  return out; 
}


// [[Rcpp::export]]
NumericVector compute_surv_m_t(NumericVector surv_prev,
                               IntegerVector first,
                               int t,
                               double phi){
  
  int n = first.size();
  NumericVector out(n, 1.0);
  
  for(int i = 0; i < n; i++){
    if(first[i] < t){
      if(surv_prev[i] == 1){
        out[i] = rbinom(1,1,phi)[0];
      } else {
        out[i] = 0;
      }
     
    }
  }
  
  return out; 
}


// [[Rcpp::export]]
NumericVector compute_surv_f_t(NumericVector surv_prev,
                               NumericVector surv_partner,
                               NumericVector Phi_cond_t,
                               double phi, 
                               int t,
                               IntegerVector first_f,
                               IntegerVector first_m,
                               IntegerVector pairs){
  
  int nf = first_f.size();
  int nm = first_m.size();
  NumericVector out(nf, 1.0);
  
  for(int i = 0; i < nf; i++){
    if(first_f[i] < t){
      
      int j = pairs[i];
      
      if((j == nm)){
        out[i] = rbinom(1,1,phi * surv_prev[i])[0];
      } else {
        
        if(surv_partner[j] == 0.0){
          out[i] = rbinom(1,1,Phi_cond_t[0] * surv_prev[i])[0];
        } else {
          out[i] = rbinom(1,1,Phi_cond_t[1] * surv_prev[i])[0]; 
        } 
      }
    } 
  }
  
  return out; 
}

double inv_logit(double eta){
  return 1/(1 + exp(-eta));
}


NumericVector compute_repartner_t(NumericVector mating_f,
                                  NumericVector mating_m,
                                  IntegerMatrix pairs,
                                  NumericVector beta,
                                  int t){
  
  int nf = mating_f.size();
  int nm = mating_m.size();
  NumericVector out(nf);
  
  for(int i = 0; i < nf; i++){
    
    int j = pairs(i,t-1);
    
    if(j != nm){
      
      double hij = 1;
      
      for(int h = 0; h < (t-1); h++){
        if(j == pairs(i,h)){
          hij += 1; 
        }
      }
      
      double eta = beta[0]  + beta[1] * hij;
      double prob = inv_logit(eta) * mating_f[i] * mating_m[j];
      out[i] = rbinom(1,1,prob)[0];
    }
  }
  return out;
}

NumericMatrix compute_recapture_m(NumericMatrix surv,
                                  IntegerVector first, 
                                  int k,
                                  NumericVector p){
  
  int nm = first.size();
  NumericMatrix xm(nm, k);
  
  for(int i = 0; i < nm; i++){
    xm(i,first[i]) = 1; 
    for(int t = (first[i]+1); t < k; t++){
      xm(i,t) = rbinom(1,1,p[t] * surv(i, t))[0];
    }
  }
  
  return xm;
}

NumericMatrix compute_recapture_f(NumericMatrix surv,
                                  NumericMatrix capture_partner,
                                  NumericMatrix P_cond,
                                  NumericVector p, 
                                  int k,
                                  IntegerVector first_f,
                                  IntegerVector first_m,
                                  IntegerMatrix pairs){
  
  
  int nf = first_f.size();
  int nm = first_m.size();
  NumericMatrix xf(nf, k);
  
  
  for(int i = 0; i < nf; i++){
    xf(i,first_f[i]) = 1; 
    for(int t = (first_f[i]+1); t < k; t++){
      int j = pairs(i,t);
      bool single = (j == nm);
      if(single){
        xf(i,t) = rbinom(1,1,p[t] * surv(i,t))[0];
      } else {
        xf(i,t) = rbinom(1,1,P_cond(capture_partner(j,t), t) * surv[i])[0];
      }
    }
  }
  
  return xf;
}

NumericMatrix compute_hidden_survival(NumericMatrix x,
                                      IntegerVector first,
                                      NumericMatrix s, 
                                      int k){
  int n = first.size();
  NumericMatrix a(n,k);
  NumericVector last(n);
  
  for(int i = 0; i < n; i++){
    for(int t= 0; t < k; t++){
      if(x(i,t) == 1){
        last[i] = t;
      }
    }
    
    for(int t= first[i]; t < (1+last[i]); t++){
      a(i,t) = 1;
    }
    
    if(last[i] !=(k-1)){
      for(int t= (1+last[i]); t < k; t++){
        a(i,t) = NumericVector::get_na();
      }
    }
    
  }
  
  return a;
}

NumericMatrix compute_hidden_pairs(IntegerMatrix pairs_f,
                                   NumericMatrix xf,
                                   NumericMatrix xm,
                                   NumericMatrix af,
                                   NumericMatrix am,
                                   NumericMatrix repartner,
                                   NumericMatrix mating_f,
                                   NumericMatrix mating_m,
                                   IntegerVector first_f,
                                   IntegerVector first_m,
                                   int k,
                                   bool imputed_pairs){
  
  int nf = first_f.size();
  int nm = first_m.size();
  NumericMatrix apairs_f(nf, k);
  
  NumericVector last_f(nf);
  
  for(int i = 0; i < nf; i++){
    for(int t= 0; t < k; t++){
      if(xf(i,t) == 1){
        last_f[i] = t;
      }
    }
  }
  
  for(int i = 0; i < nf; i++){
    for(int t = 0; t < k; t++){
      if(first_f[i] > t){
        apairs_f(i,t) = nm;
      } else {
       int j = pairs_f(i,t);
        if(j != nm){
          if((xm(j,t) == 1) & (xf(i,t) == 1)){
            apairs_f(i,t) = j;
          } else {
            apairs_f(i,t) = NumericMatrix::get_na();
          }
        } else {
          if(xf(i,t) == 1){
            apairs_f(i,t) = nm;
          } else {
            apairs_f(i,t) = NumericMatrix::get_na();
          }
        }
      }
    }
  }
  
  //Imputing NAs between pairs
  if(imputed_pairs){
    for(int i = 0; i < nf; i++){
      
      for(int t = (first_f[i]+1); t < (last_f[i]+1); t++){
        int j_prev = apairs_f(i,t-1);
        int j = apairs_f(i,t);
        
        // If partner is unknown
        if(IntegerVector::is_na(j)){
          if(!IntegerVector::is_na(j_prev)){
            if(j_prev != nm){
              if(!(NumericVector::is_na(am(j_prev,t))) & !(NumericVector::is_na(af(i,t)))){
                apairs_f(i,t) = j_prev;
              }
            }
          }
        }
      }
    }
  }
  
  
  // Update output to match R indexing
    for(int i = 0; i < nf; i++){
      
      for(int t = 0; t < k; t++){
        int j = apairs_f(i,t);
        
        // If not unknown add 1
        if(!IntegerVector::is_na(j)){
          apairs_f(i,t) = apairs_f(i,t) + 1;
        }
      }
  }
  
  return apairs_f;

}


// [[Rcpp::export]]
List simulate_ps_data(int n,
                      int k,
                      double propfemale,
                      NumericVector delta,
                      NumericVector phif,
                      NumericVector phim,
                      NumericVector gam,
                      NumericVector pf,
                      NumericVector pm,
                      NumericVector rho,
                      NumericVector beta,
                      NumericVector tau,
                      bool imputed_pairs){

  //Probability Matrices
  NumericMatrix Phi      = compute_joint_prob(phif, phim, gam);
  NumericMatrix P        = compute_joint_prob(pf,   pm,   rho);
  NumericMatrix Phi_Cond = compute_cond_prob(Phi,   phif, phim);
  NumericMatrix P_Cond   = compute_cond_prob(P,     pf,   pm);

  // COnstruct Data Objects

  // Males and Females
  int nf = floor(n* propfemale);
  int nm = n - nf;

  // Mark-Recapture Objects

  //Female capture/recruitment/mating/survival
  IntegerVector first_capture_f = initialize_first_capture(nf, k , tau)-1;
  NumericMatrix mating_f        = initialize_mating(nf, k, first_capture_f, delta);
  NumericMatrix sf              = initialize_survival(nf, k, first_capture_f);
  NumericMatrix repartner(nf, k);

  //Male capture/recruitment/mating/survival
  IntegerVector first_capture_m = initialize_first_capture(nm, k , tau)-1;
  NumericMatrix mating_m        = initialize_mating(nm, k, first_capture_m, delta);
  NumericMatrix sm              = initialize_survival(nm, k, first_capture_m);

  // Start Condition for looping over time
  int t0f = minIntCppV(first_capture_f);
  int t0m = minIntCppV(first_capture_m);
  int t0  = minIntCpp(t0f, t0m);

  // Initialize Pairs Matrix
  IntegerMatrix pairs_f         = initialize_pairs(nf, nm, k, t0, mating_f, mating_m);

  // // Generate MR data iteratively on time
  for(int t = t0+1; t < k; t++){
      sm( _ , t)        = compute_surv_m_t(sm( _ , t-1), first_capture_m, t, phim[t-1]);
      sf( _ , t)        = compute_surv_f_t(sf( _ , t-1),
                                           sm( _ , t),
                                           Phi_Cond( _ , t-1),
                                           phif[t-1],
                                           t,
                                           first_capture_f,
                                           first_capture_m,
                                           pairs_f( _ ,t-1));

      mating_m( _ , t)  = compute_mating_t(sm( _ , t), first_capture_m, t, delta[t]);
      mating_f( _ , t)  = compute_mating_t(sf( _ , t), first_capture_f, t, delta[t]);
      repartner( _ , t) = compute_repartner_t(mating_f(_, t),
                                              mating_m(_, t),
                                              pairs_f,
                                              beta,
                                              t);
      pairs_f( _ , t)   = compute_pairs_t(repartner( _ , t),
                                          pairs_f(_, t-1),
                                          mating_f(_, t),
                                          mating_m(_, t),
                                          nf,
                                          nm);
  }

  // // Compute Recapture Histories
  NumericMatrix xm = compute_recapture_m(sm, first_capture_m, k, pm);
  NumericMatrix xf = compute_recapture_f(sf, xm, P_Cond, pf, k, first_capture_f, first_capture_m, pairs_f);

  // // Compute Hidden Survival
  NumericMatrix af = compute_hidden_survival(xf, first_capture_f, sf, k);
  NumericMatrix am = compute_hidden_survival(xm, first_capture_m, sm, k);

  // // Compute Hidden Pairs
  NumericMatrix apairs_f = compute_hidden_pairs(pairs_f,
                                                xf,
                                                xm,
                                                af,
                                                am,
                                                repartner,
                                                mating_f,
                                                mating_m,
                                                first_capture_f,
                                                first_capture_m,
                                                k,
                                                imputed_pairs);

  List out(14);
  CharacterVector list_names(14);
  
  list_names[0] = "n";
  list_names[1] = "k";
  list_names[2] = "nf";
  list_names[3] = "nm";
  list_names[4] = "first_capture_f";
  list_names[5] = "first_capture_m";
  list_names[6] = "pairs_f";
  list_names[7] = "sf";
  list_names[8] = "sm";
  list_names[9] = "apairs_f";
  list_names[10] = "af";
  list_names[11] = "am";
  list_names[12] = "recap_f";
  list_names[13] = "recap_m";
  
  out[0]  = n;
  out[1]  = k;
  out[2]  = nf;
  out[3]  = nm;
  out[4]  = first_capture_f + 1;
  out[5]  = first_capture_m + 1;
  out[6]  = pairs_f;
  out[7]  = sf;
  out[8]  = sm;
  out[9]  = apairs_f;
  out[10] = af;
  out[11] = am;
  out[12] = xf;
  out[13] = xm;
  
  out.attr("names") = list_names;
  
  return out;
}