var x{1..6} integer >= 0 <= 3;

maximize o: -x[1] + 0*x[2] + 0*x[3] + 2*x[4] - 5*x[5] - 3*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: abs(x[5] - x[4]) <= 0;
s.t. c3: (if x[1] > 0 then x[1] else x[2]) < 2;
