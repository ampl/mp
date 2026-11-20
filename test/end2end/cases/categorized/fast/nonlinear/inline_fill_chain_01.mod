# -------------------------------------------------------------
# inline_fill_chain_01.mod
# Chain y[k] = y[k-1] + count(...): each COUNT is linearized during
# conversion; inlining the shared chain links into each other would grow
# quadratically. cvt:pre:linfill (default 10) keeps the links as variables
# once the fill bound exceeds the limit.
# -------------------------------------------------------------
param N := 20; param M := 3;
param a{k in 1..N, j in 1..M} := ((7*k + 3*j) mod 9)/2 - 2;
var x{1..N, 1..M} >= -3, <= 3;
var y{k in 0..N} = if k == 0 then 0 else y[k-1] + count{j in 1..M} (x[k,j] >= a[k,j]);
minimize o: sum{k in 1..N} exp(y[k]/(N*M)) - sum{k in 1..N, j in 1..M} 0.01*x[k,j];
