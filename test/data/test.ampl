var x{i in 1..3} >= 10 + i, <= 20 + i;
var y{i in 1..2} >= 13 + i, <= 23 + i, integer;

minimize o1{i in 1..16}:
  (30 + i) * x[1] + sin(sum{j in 1..3} x[j] + sum{j in 1..2} y[j]);
minimize o2: cos(x[1]);
maximize o3{i in 1..2}: (50 + i) * y[1];

s.t. c1{i in 1..10}:
  100 + i <= (60 + i) * x[1] +
             log(sum{j in 1..3} x[j] + sum{j in 1..2} y[j]) <= 200 + i;
s.t. c2: exp(x[1]) == 0;
s.t. c3{i in 1..2}: 111 + i <= (80 + i) * x[3] <= 211 + i;
s.t. c4{1..6}: x[1] != 1;
s.t. c5: x[1] != 1 && x[2] != 2;

option presolve 0;
