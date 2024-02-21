var x {1..7} >= 0;
var y >= 0;

maximize obj: sum {j in 1..7} -x[j]^2 - y^2;

subj to conQ: sum {j in 1..7} x[j]^2 <= y^2;
subj to conL: sum {j in 1..7} x[j] = 75;
