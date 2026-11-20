var x{1..5} integer >= 0 <= 2;

maximize o: -3*x[1] + 2*x[2] + x[3] - 2*x[4] - x[5];

s.t. c1: (x[1] > 1) or (x[2] != 0) or (x[2] >= 1);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: ((if x[1] == 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c4: (if x[3] != 0 then x[2] else x[4]) < 1;
