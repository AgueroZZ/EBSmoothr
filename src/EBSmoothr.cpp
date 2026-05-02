#include <TMB.hpp>

template<class Type>
Type log_matern_pc_prior(Type log_range,
                         Type log_sigma,
                         vector<Type> pc_prior,
                         int d)
{
  Type d_type = Type(d);
  Type range = exp(log_range);
  Type sigma = exp(log_sigma);

  Type rho0 = pc_prior(0);
  Type alpha_r = pc_prior(1);
  Type sigma0 = pc_prior(2);
  Type alpha_s = pc_prior(3);

  Type R = -log(alpha_r) * pow(rho0, d_type / Type(2.0));
  Type S = -log(alpha_s) / sigma0;

  Type log_prior_range =
    log(d_type * R / Type(2.0)) -
    (d_type / Type(2.0)) * log_range -
    R * pow(range, -d_type / Type(2.0));
  Type log_prior_sigma = log(S) + log_sigma - S * sigma;

  return log_prior_range + log_prior_sigma;
}

template<class Type>
Type log_noise_pc_prior(Type log_noise,
                        vector<Type> pc_prior)
{
  Type noise0 = pc_prior(4);
  Type alpha_noise = pc_prior(5);
  Type noise = exp(log_noise);
  Type S = -log(alpha_noise) / noise0;

  return log(S) + log_noise - S * noise;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(model_id);          // 0: L-GP, 1: Matern

  if (model_id == 0) {
    // --------------------
    // data
    // --------------------
    DATA_VECTOR(x);                // observed
    DATA_VECTOR(s);                // known sd (length n)
    DATA_SPARSE_MATRIX(X);         // global design (n x pX)
    DATA_SPARSE_MATRIX(B);         // local design  (n x pB)
    DATA_SPARSE_MATRIX(P);         // penalty precision (pB x pB), SPD
    DATA_SCALAR(logPdet);          // log det(P), precomputed
    DATA_SCALAR(betaprec);         // beta precision; <=0 => diffuse/no prior
    DATA_INTEGER(link_id);         // 0: identity, 1: exp-link (mu = exp(eta))
    DATA_INTEGER(learn_noise);     // 0: use known s, 1: learn one common noise SD

    // --------------------
    // parameters
    // --------------------
    PARAMETER(theta);              // scalar; controls U prior scale
    PARAMETER_VECTOR(U);           // length pB
    PARAMETER_VECTOR(beta);        // length pX
    PARAMETER(log_noise);          // scalar log common noise SD when learn_noise = 1

    // --------------------
    // linear predictor
    // --------------------
    vector<Type> eta = X * beta + B * U;

    // apply link
    vector<Type> mu = eta;
    if (link_id == 1) {
      mu = exp(eta);
    }

    vector<Type> noise_sd = s;
    if (learn_noise == 1) {
      noise_sd.setConstant(exp(log_noise));
    }

    // --------------------
    // likelihood: x ~ N(mu, s^2)
    // --------------------
    Type ll = sum(dnorm(x, mu, noise_sd, true));

    // --------------------
    // prior on U: Gaussian with precision exp(theta) * P
    // --------------------
    Type lp = Type(0);

    int pB = P.cols();
    if (pB > 0) {
      vector<Type> PU = P * U;
      Type UPU = (U * PU).sum();

      lp += -Type(0.5) * exp(theta) * UPU;
      lp += Type(0.5) * (Type(pB) * theta + logPdet) -
        Type(0.5) * Type(pB) * log(Type(2.0) * M_PI);
    }

    // --------------------
    // prior on beta only if betaprec > 0 (proper)
    // --------------------
    int pX = X.cols();
    if (betaprec > Type(0) && pX > 0) {
      Type bb = (beta * beta).sum();
      lp += -Type(0.5) * betaprec * bb;
      lp += Type(0.5) * Type(pX) * log(betaprec) -
        Type(0.5) * Type(pX) * log(Type(2.0) * M_PI);
    }

    return -(ll + lp);
  }

  if (model_id == 1) {
    // --------------------
    // data
    // --------------------
    DATA_VECTOR(x);                // observed
    DATA_VECTOR(s);                // known sd (length n)
    DATA_SPARSE_MATRIX(A);         // observation projector (n x n_spde)
    DATA_SPARSE_MATRIX(M0);        // SPDE precision basis
    DATA_SPARSE_MATRIX(M1);        // SPDE precision basis
    DATA_SPARSE_MATRIX(M2);        // SPDE precision basis
    DATA_SCALAR(betaprec);         // beta precision; <=0 => no proper prior
    DATA_SCALAR(matern_alpha);     // v1 expects alpha = 2
    DATA_INTEGER(matern_d);        // spatial dimension
    DATA_INTEGER(link_id);         // 0: identity, 1: exp-link (mu = exp(eta))
    DATA_INTEGER(learn_noise);     // 0: use known s, 1: learn one common noise SD
    DATA_INTEGER(use_pc_prior);    // 0: no PC prior, 1: include PC prior
    DATA_INTEGER(use_pc_noise_prior); // 0: no noise PC prior, 1: include noise PC prior
    DATA_VECTOR(pc_prior);         // range anchor/alpha, sigma anchor/alpha, optional noise anchor/alpha

    // --------------------
    // parameters
    // --------------------
    PARAMETER(log_range);
    PARAMETER(log_sigma);
    PARAMETER_VECTOR(w);           // latent SPDE weights
    PARAMETER_VECTOR(beta);        // scalar intercept, kept as vector for mapping
    PARAMETER(log_noise);          // scalar log common noise SD when learn_noise = 1

    Type nll = Type(0);

    // --------------------
    // Matern SPDE precision for alpha = 2:
    // Q = tau^2 * (kappa^4 M0 + 2 kappa^2 M1 + M2)
    // --------------------
    Type d_type = Type(matern_d);
    Type nu = matern_alpha - d_type / Type(2.0);
    Type log_kappa = Type(0.5) * log(Type(8.0) * nu) - log_range;
    Type kappa2 = exp(Type(2.0) * log_kappa);
    Type kappa4 = kappa2 * kappa2;
    Type log_tau =
      Type(0.5) * (
        lgamma(nu) -
        lgamma(matern_alpha) -
        (d_type / Type(2.0)) * log(Type(4.0) * M_PI) -
        Type(2.0) * nu * log_kappa -
        Type(2.0) * log_sigma
      );
    Type tau2 = exp(Type(2.0) * log_tau);

    Eigen::SparseMatrix<Type> Q = M0 * kappa4;
    Q += M1 * (Type(2.0) * kappa2);
    Q += M2;
    Q *= tau2;

    // --------------------
    // likelihood: x ~ N(mu, s^2)
    // --------------------
    vector<Type> eta = A * w;
    for (int i = 0; i < eta.size(); i++) {
      eta(i) += beta(0);
    }

    vector<Type> mu = eta;
    if (link_id == 1) {
      mu = exp(eta);
    }
    vector<Type> noise_sd = s;
    if (learn_noise == 1) {
      noise_sd.setConstant(exp(log_noise));
    }
    nll -= sum(dnorm(x, mu, noise_sd, true));

    // --------------------
    // latent and beta priors
    // --------------------
    nll += density::GMRF(Q)(w);

    if (betaprec > Type(0)) {
      nll -= dnorm(beta(0), Type(0), Type(1.0) / sqrt(betaprec), true);
    }

    if (use_pc_prior == 1) {
      nll -= log_matern_pc_prior(log_range, log_sigma, pc_prior, matern_d);
    }
    if (learn_noise == 1 && use_pc_noise_prior == 1) {
      nll -= log_noise_pc_prior(log_noise, pc_prior);
    }

    return nll;
  }

  Rf_error("Unsupported model_id in the EBSmoothr TMB objective.");
  return Type(0);
}
