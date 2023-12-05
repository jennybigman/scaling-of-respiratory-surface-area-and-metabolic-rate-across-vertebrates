data {
int N; // sample size
vector[N] RSA; // response (respiratory surface area)
vector[N] MRSA; // predictor (body mass)
matrix[N, N] d_mat; // sigma matrix
matrix[N, N] vcov_mat; // vcov matrix
}

parameters {
real a; // intercept
real bMass; // mass-slope
real <lower=0> sigma; // error
real <lower=0,upper=1> lambda; // phylogenetic signal
}

transformed parameters {

matrix[N, N] sigma_mat;
matrix[N, N] sigma_total;

vector[N] mu_RSA;
vector[N] resid;

sigma_mat = (1-lambda)*d_mat + lambda*vcov_mat;
sigma_total = sigma*sigma_mat;

for(i in 1:N) {
mu_RSA[i] = a + bMass * MRSA[i];
resid[i] = RSA[i] - mu_RSA[i];
}}


model {
lambda  ~ uniform(0,1);
a       ~ student_t(3, 0, 10);
bMass   ~ student_t(3, 0, 10);
sigma   ~ cauchy(0, 10);

RSA ~ multi_normal(mu_RSA, sigma_total);
}

generated quantities {
real log_lik;
log_lik = multi_normal_lpdf(RSA | mu_RSA, sigma_total);
}
