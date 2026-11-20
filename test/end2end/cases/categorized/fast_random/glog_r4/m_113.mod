var x{1..3} integer >= 0 <= 3;

minimize o: -x[1] - 2*x[2] - 5*x[3];

s.t. c1: ((if x[2] == 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: abs(x[2] - x[3]) <= 2;
