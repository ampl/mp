var x{1..6} integer >= 0 <= 2;

maximize o: 4*x[1] - 3*x[2] + 3*x[3] - 4*x[4] + 0*x[5] - 5*x[6];

s.t. c1: abs(x[2] - x[2]) < 0;
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) > 1;
s.t. c4: (if x[5] > 2 then x[6] else x[6]) >= 0;
