# -------------------------------------------------------------
# inline_quad_bin_01.mod
# Products of binaries inside nested linear expressions and indicators
# -------------------------------------------------------------
var u{1..4} binary; var b binary;
minimize o: -sum{i in 1..4} u[i] + 2*b;
s.t. c1: b ==> u[1]*u[2] + u[3]*u[4] + (if not (u[1] + u[3] >= 1) then 1) <= 1;
s.t. c2: (not b) ==> u[1] + u[2] <= 1;
s.t. c3: u[2]*u[3] + (if u[4] >= 1 then 1) <= 1;
