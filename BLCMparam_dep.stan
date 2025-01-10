data {
  int<lower=1> P;               // Number of populations
  int<lower=1> n[P];            // Sample size for each population
  int<lower=1> N;               // Total number of observations
  int<lower=1, upper=P> Pop[N]; // Population indicator
  int<lower=0, upper=1> Test1[N]; // Binary outcome of test 1
  int<lower=0, upper=1> Test2[N]; // Binary outcome of test 2
}

parameters {
  real<lower=0, upper=1> pi[P];    // Prevalence for each population
  real<lower=0, upper=1> Se1;      // Sensitivity of test 1
  real<lower=0, upper=1> Sp1;      // Specificity of test 1
  real<lower=0, upper=1> Se2;      // Sensitivity of test 2
  real<lower=0, upper=1> Sp2;      // Specificity of test 2
  real<lower=-1, upper=1> cdp;     // Correlation for positive tests
  real<lower=-1, upper=1> cdn;     // Correlation for negative tests
}

model {
  // Priors for sensitivity, specificity, dependences, and prevalence
  Se1 ~ beta(3, 2);  // Sensitivity for test 1
  Sp1 ~ beta(3, 2);  // Specificity for test 1
  Se2 ~ beta(3, 2);  // Sensitivity for test 2
  Sp2 ~ beta(3, 2);  // Specificity for test 2

  for (p in 1:P) {
    pi[p] ~ beta(2, 3);  // Prior for prevalence in each population
  }

  // Likelihood
  for (i in 1:N) {
    real prob11 = fmax(pi[Pop[i]] * (Se1 * Se2 + cdp), 1e-10) + fmax((1 - pi[Pop[i]]) * ((1 - Sp1) * (1 - Sp2) + cdn), 1e-10);
    real prob10 = fmax(pi[Pop[i]] * (Se1 * (1 - Se2) - cdp), 1e-10) + fmax((1 - pi[Pop[i]]) * ((1 - Sp1) * Sp2 - cdn), 1e-10);
    real prob01 = fmax(pi[Pop[i]] * ((1 - Se1) * Se2 - cdp), 1e-10) + fmax((1 - pi[Pop[i]]) * (Sp1 * (1 - Sp2) - cdn), 1e-10);
    real prob00 = fmax(pi[Pop[i]] * ((1 - Se1) * (1 - Se2) + cdp), 1e-10) + fmax((1 - pi[Pop[i]]) * (Sp1 * Sp2 + cdn), 1e-10);

    if (Test1[i] == 1 && Test2[i] == 1) {
      target += log(prob11);  // Log-likelihood for (1, 1)
    } else if (Test1[i] == 1 && Test2[i] == 0) {
      target += log(prob10);  // Log-likelihood for (1, 0)
    } else if (Test1[i] == 0 && Test2[i] == 1) {
      target += log(prob01);  // Log-likelihood for (0, 1)
    } else if (Test1[i] == 0 && Test2[i] == 0) {
      target += log(prob00);  // Log-likelihood for (0, 0)
    }
  }
}
