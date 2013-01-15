var x{i in 1..2} >= 13 + i, <= 23 + i, integer;
var y{i in 1..3} >= 10 + i, <= 20 + i;

minimize o1{1..17}: sin(sum{i in 1..2} x[i] + sum{i in 1..3} y[i]);
maximize o2{i in 1..2}: x[i] + y[i];

s.t. c1{i in 1..11}: 100 + i <= sin(sum{j in 1..2} x[j] + sum{j in 1..3} y[j]) <= 200 + i;
s.t. c2{i in 1..2}: 111 + i <= x[i] + y[i] <= 211 + i;
s.t. c3{1..7}: x[1] != 1;

option presolve 0;
