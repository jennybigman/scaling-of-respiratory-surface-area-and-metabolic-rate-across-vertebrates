// mod 1

data {
int N; // sample size
vector[N] MR; // response (metabolic rate)
vector[N] MMR; // predictor (body mass)
vector[N] InTemp; // predictor (temperature)
matrix[N, N] d_mat; // sigma matrix
matrix[N, N] vcov_mat; // vcov matrix
}

parameters {
real a; // intercept 
real bMass; // slope of body mass
real bTemp; // slope of temperature
real <lower=0> sigma; // error
real <lower=0,upper=1>lambda; // phylogenetic signal
}

transformed parameters {

matrix[N, N] sigma_mat;
matrix[N, N] sigma_total;

vector[N] mu_MR;

sigma_mat = (1-lambda)*d_mat + lambda*vcov_mat;
sigma_total = sigma_mat * sigma;

for(i in 1:N) {
mu_MR[i] = a + bMass * MMR[i] + bTemp * InTemp[i];

}
}


model {

lambda  ~ uniform(0,1);
a       ~ student_t(3, -1, 10);
bMass   ~ student_t(3, 0, 10);
bTemp   ~ student_t(3, 0, 10);
sigma   ~ cauchy(0, 10);

MR ~ multi_normal(mu_MR, sigma_total);

}

generated quantities {

real log_lik;

log_lik =
multi_normal_lpdf(MR | mu_MR, sigma_total);

}
