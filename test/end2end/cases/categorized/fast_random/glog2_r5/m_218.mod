var x{1..6} integer >= 0 <= 1;

maximize o: 3*x[1] + 6*x[2] + x[3] + 2*x[4] - 3*x[5] + 3*x[6];

s.t. c1: (x[5] != 1) and (x[5] >= 0) ==> (x[2] >= 1);
s.t. c2: abs(x[6] - x[6]) >= 1;
