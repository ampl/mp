var x{1..4} integer >= 0 <= 2;

maximize o: -4*x[1] - 2*x[2] + 5*x[3] + x[4];

s.t. c1: (if x[3] <= 1 then x[3] else x[4]) >= 1;
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: not ((x[4] < 2) and (x[2] != 1));
s.t. c4: ((if x[1] >= 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c5: x[3] == 1;
