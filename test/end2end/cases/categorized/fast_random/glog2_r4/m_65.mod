var x{1..5} integer >= 0 <= 4;

maximize o: 0*x[1] + 3*x[2] - 4*x[3] + x[4] + 2*x[5];

s.t. c1: ((if x[1] > 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
s.t. c2: max(x[1], x[1]) - min(x[3], x[4]) >= 4;
s.t. c3: max(x[3], x[2]) - min(x[3], x[3]) >= 0;
s.t. c4: ((x[2] > 0) ==> (x[4] > 2)) and ((x[4] > 2) ==> (x[5] >= 0));
s.t. c5: forall {i in 1..4} x[i] != x[i+1];
s.t. c6: abs(x[5] - x[2]) != 1;
s.t. c7: not ((x[4] > 0) and (x[5] > 0));
