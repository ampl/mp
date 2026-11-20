var x{1..5} integer >= 0 <= 3;

maximize o: -3*x[1] - 2*x[2] + 2*x[3] + 0*x[4] + 5*x[5];

s.t. c1: (if x[5] < 0 then x[5] else x[4]) >= 0;
s.t. c2: atleast 2 {i in 1..5} (x[i] >= 3);
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
