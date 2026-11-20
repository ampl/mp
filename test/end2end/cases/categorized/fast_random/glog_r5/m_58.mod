var x{1..6} integer >= 0 <= 4;

maximize o: -2*x[1] + 3*x[2] - 2*x[3] - 5*x[4] + 4*x[5] - 5*x[6];

s.t. c1: (if x[3] > 1 then x[1] else x[6]) >= 1;
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: alldiff(x[6], x[1], x[5], x[2]);
