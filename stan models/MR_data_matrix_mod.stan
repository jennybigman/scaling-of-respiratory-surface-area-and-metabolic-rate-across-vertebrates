      data {
        int N; // sample size
        vector[N] MR; // response (metabolic rate)
        int K; // number of predictors
        matrix[N, K] data_matrix; // matrix of predictors
        matrix[N, N] d_mat; // sigma matrix
        matrix[N, N] vcov_mat; // vcov matrix
      }
      
      parameters {
        vector[K] beta; // coefficients
        real<lower=0> sigma; // error
        real<lower=0,upper=1> lambda; // phylogenetic signal
      }
      
      
      transformed parameters {
        
        matrix[N, N] sigma_mat;
        matrix[N, N] sigma_total;
        
        vector[N] mu_MR;
        

        sigma_mat = (1-lambda)*d_mat + lambda*vcov_mat;
        sigma_total = sigma*sigma_mat;
        
        mu_MR = data_matrix * beta;

      }
      
      
      model {
        
        
        lambda  ~ uniform(0,1);
        beta    ~ student_t(3, 0, 10);
        sigma   ~ cauchy(0, 10);
        
        
        MR ~ multi_normal(mu_MR, sigma_total);
        
      }
      
      
      generated quantities {
        
        real log_lik;
        
        log_lik = multi_normal_lpdf(MR | mu_MR, sigma_total); 
        
      }
