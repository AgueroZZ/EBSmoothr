#include <TMB.hpp>
//#include <fenv.h>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // --------------------
  // data
  // --------------------
  DATA_VECTOR(x);                  // observed
  DATA_VECTOR(s);                  // known sd (length n)
  DATA_SPARSE_MATRIX(X);           // global design (n x pX)
  DATA_SPARSE_MATRIX(B);           // local design  (n x pB)
  DATA_SPARSE_MATRIX(P);           // penalty precision (pB x pB), SPD
  DATA_SCALAR(logPdet);            // log det(P), precomputed
  DATA_SCALAR(betaprec);           // beta precision; <=0 => diffuse/no prior
  DATA_INTEGER(link_id);           // 0: identity, 1: exp-link (mu = exp(eta))

  // --------------------
  // parameters
  // --------------------
  PARAMETER(theta);                // scalar; controls U prior scale
  PARAMETER_VECTOR(U);             // length pB
  PARAMETER_VECTOR(beta);          // length pX

  // --------------------
  // linear predictor
  // --------------------
  vector<Type> eta = X * beta + B * U;   // n-vector

  // apply link
  vector<Type> mu = eta;
  if (link_id == 1) {
    mu = exp(eta);
  }

  // --------------------
  // likelihood: x ~ N(mu, s^2)
  // --------------------
  Type ll = sum(dnorm(x, mu, s, true));

  // --------------------
  // prior on U: Gaussian with precision exp(theta) * P
  //   U' (exp(theta) P) U / 2  + normalizing constant
  // --------------------
  Type lp = Type(0);

  int pB = P.cols();
  if (pB > 0) {
    vector<Type> PU = P * U;
    Type UPU = (U * PU).sum();

    // quadratic form
    lp += -Type(0.5) * exp(theta) * UPU;

    // log normalizing constant: + 0.5 * logdet(exp(theta) P) - 0.5*pB*log(2pi)
    // logdet(exp(theta) P) = pB*theta + logdet(P)
    lp += Type(0.5) * (Type(pB) * theta + logPdet) - Type(0.5) * Type(pB) * log(Type(2.0) * M_PI);
  }

  // --------------------
  // prior on beta only if betaprec > 0 (proper)
  // beta ~ N(0, (betaprec I)^{-1})
  // --------------------
  int pX = X.cols();
  if (betaprec > Type(0) && pX > 0) {
    Type bb = (beta * beta).sum();
    lp += -Type(0.5) * betaprec * bb;
    lp += Type(0.5) * Type(pX) * log(betaprec) - Type(0.5) * Type(pX) * log(Type(2.0) * M_PI);
  }

  // negative log posterior
  return -(ll + lp);
}
