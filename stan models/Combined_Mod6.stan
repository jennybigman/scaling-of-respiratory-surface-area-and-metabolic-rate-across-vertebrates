
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

parameters {
real aRSA;
real bMassRSA;
real<lower=0> sigma_RSA; // error for RSA model

real aMR;
real bMassMR;
real bResid;
real bIntMassTS;
real bIntMassResid;
real bTemp;
real bTherm;
real<lower=0> sigma_MR; // error for MR model

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

sigma_mat_MR = (1-lambda_MR)*d_mat + lambda_MR*vcov_mat;
sigma_total_MR = sigma_MR*sigma_mat_MR;

sigma_mat_RSA = (1-lambda_RSA)*d_mat + lambda_RSA*vcov_mat;
sigma_total_RSA = sigma_RSA*sigma_mat_RSA;


for(i in 1:N) {
mu_RSA[i] = aRSA + bMassRSA * MRSA[i];
resid[i] = RSA[i] - mu_RSA[i];
mu_MR[i] = aMR + bResid * resid[i] + bMassMR * MMR[i] + bTherm * ThermoStrat[i] + bTemp * Temp[i] + bIntMassResid * MMR[i] * resid[i] + bIntMassTS * MMR[i] * ThermoStrat[i];
}}

model {


aRSA        ~ student_t(3, 0, 10);
bMassRSA    ~ student_t(3, 0, 10);
aMR         ~ student_t(3, 0, 10);
bMassMR     ~ student_t(3, 0, 10);
bTemp       ~ student_t(3, 0, 10);
bTherm      ~ student_t(3, 0, 10);
bIntMassResid ~ student_t(3, 0, 10);
bIntMassTS ~ student_t(3, 0, 10);
sigma_MR    ~ cauchy(0, 10);
sigma_RSA   ~ cauchy(0, 10);
lambda_MR   ~ uniform(0, 1);
lambda_RSA  ~ uniform(0, 1);


RSA ~ multi_normal(mu_RSA, sigma_total_RSA);
MR ~ multi_normal(mu_MR, sigma_total_MR);
}

generated quantities {

real log_lik;

log_lik =
multi_normal_lpdf(RSA | mu_RSA, sigma_total_RSA) +
multi_normal_lpdf(MR | mu_MR, sigma_total_MR);
}
