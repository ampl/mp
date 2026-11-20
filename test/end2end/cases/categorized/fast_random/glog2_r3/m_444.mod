var x{1..4} integer >= -3 <= 2;

maximize o: -2*x[1] - 2*x[2] + 3*x[3] - 4*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: (if x[4] == 0 then x[2] else x[2]) < 2;
s.t. c3: (x[1] < 0) or (x[1] == -3) or (x[1] != 1);
