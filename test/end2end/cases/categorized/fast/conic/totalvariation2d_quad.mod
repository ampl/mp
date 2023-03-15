## 2D Total Variation, using quadratic formulation.
## This case study is based mainly on the paper by Goldfarb and Yin [GY05]
## https://docs.mosek.com/latest/pythonfusion/case-studies-total-variation.html

param n>1;
param m>1;

param f {i in 1..n, j in 1..m} default        ## signal + noise
    1.0*(i+j-2)/(n+m) + Normal(0.0, 0.08);
param sigmaVal default 0.0004;
param sigma >=0.0 default sigmaVal*n*m;

var u {1..n, 1..m} >=0.0, <=1.0;
var delta_uf {1..n, 1..m};
var t {1..n-1, 1..m-1} >=0.0;
var delta_plus {1..n-1, 1..m-1, 1..2};

minimize TotalVar:
    sum {i in 1..n-1, j in 1..m-1} t[i, j];

s.t. DefDeltaPlusUp {i in 1..n-1, j in 1..m-1}:
    delta_plus[i, j, 1] == u[i+1, j] - u[i, j];

s.t. DefDeltaPlusRight {i in 1..n-1, j in 1..m-1}:
    delta_plus[i, j, 2] == u[i, j+1] - u[i, j];

s.t. DefDeltaUF {i in 1..n, j in 1..m}:
    delta_uf[i, j] == u[i, j] - f[i, j];

s.t. VariationBound {i in 1..n-1, j in 1..m-1}:
    t[i, j]^2 >= delta_plus[i, j, 1]^2 + delta_plus[i, j, 2]^2;

s.t. Total2Norm:
    -sigma^2 <= -sum {i in 1..n, j in 1..m} delta_uf[i, j]^2;

