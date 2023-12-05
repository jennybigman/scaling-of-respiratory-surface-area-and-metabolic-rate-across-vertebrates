data {
  int N; // sample size
  vector[N] RSA; // response 1 (respiratory surface area)
  vector[N] MRSA; // predictor (body mass RSA)
  vector[N] MR; // response 2 (metabolic rate)
  vector[N] MMR; // predictor (body mass MR)
  vector[N] Temp; // predictor (temperature MR)
  vector[N] ThermoStrat; // predictor (thermoregulatory strategy)
  matrix[N, N] d_mat; // sigma matrix
  matrix[N, N] vcov_mat; // vcov matrix
}

transformed data {
  vector[N] std_MMR;
  vector[N] std_Temp;
  
  vector[N] std_MRSA;
  
  std_MMR = (MMR - mean(MMR)) / sd(MMR);
  std_Temp = (Temp - mean(Temp)) / sd(Temp);
  
  std_MRSA = (MRSA - mean(MRSA)) / sd(MRSA);
}

parameters {
real aRSA;
real std_bMassRSA;
real<lower=0> std_sigma_RSA; // error for RSA model

real aMR;
real std_bMassMR;
real std_bResid;
real bIntMassTS;
real std_bTemp;
real bTherm;
real<lower=0> std_sigma_MR; // error for MR model

real<lower=0,upper=1> lambda_MR; // phylogenetic signal MR model
real<lower=0,upper=1> lambda_RSA; // phylogenetic signal RSA model
}

transformed parameters {

matrix[N, N] sigma_mat_MR;
matrix[N, N] sigma_total_MR;

matrix[N, N] sigma_mat_RSA;
matrix[N, N] sigma_total_RSA;

vector[N] mu_RSA;
vector[N] mu_MR;
vector[N] resid;
real mean_resid;
real sd_resid;
vector[N] std_resid;

sigma_mat_MR = (1-lambda_MR)*d_mat + lambda_MR*vcov_mat;
sigma_total_MR = std_sigma_MR*sigma_mat_MR;

sigma_mat_RSA = (1-lambda_RSA)*d_mat + lambda_RSA*vcov_mat;
sigma_total_RSA = std_sigma_RSA*sigma_mat_RSA;


for(i in 1:N) {
mu_RSA[i] = aRSA + std_bMassRSA * std_MRSA[i];
resid[i] = RSA[i] - mu_RSA[i];
}

mean_resid = mean(resid);
sd_resid = sd(resid);

for(i in 1:N) {
std_resid[i] = ((resid[i] - mean_resid) / sd_resid);

mu_MR[i] = aMR + std_bResid * std_resid[i] + std_bMassMR * std_MMR[i] + bTherm * ThermoStrat[i] + bIntMassTS * std_MMR[i] * ThermoStrat[i] + std_bTemp * std_Temp[i];
}}

model {

aRSA            ~ student_t(3, 0, 10);
std_bMassRSA    ~ student_t(3, 0, 10);
aMR             ~ student_t(3, 0, 10);
std_bMassMR     ~ student_t(3, 0, 10);
std_bTemp       ~ student_t(3, 0, 10);
bTherm          ~ student_t(3, 0, 10);
std_sigma_MR    ~ cauchy(0, 10);
std_sigma_RSA   ~ cauchy(0, 10);
lambda_MR       ~ uniform(0, 1);
lambda_RSA      ~ uniform(0, 1);


RSA ~ multi_normal(mu_RSA, sigma_total_RSA);
MR ~ multi_normal(mu_MR, sigma_total_MR);
}

generated quantities {

real log_lik;

log_lik =
multi_normal_lpdf(RSA | mu_RSA, sigma_total_RSA) +
multi_normal_lpdf(MR | mu_MR, sigma_total_MR);
}

