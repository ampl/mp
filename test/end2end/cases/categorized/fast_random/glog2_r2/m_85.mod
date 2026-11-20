var x{1..4} integer >= -2 <= 1;

minimize o: 5*x[1] + 5*x[2] - 4*x[3] - x[4];

s.t. c1: alldiff(x[3], x[4]);
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: ((if x[3] != -1 then 1 else 0) + (if x[4] >= -1 then 1 else 0)) >= 1;
s.t. c4: max(x[3], x[1]) - min(x[1], x[4]) <= -1;
s.t. c5: abs(x[2] - x[1]) > 1;
s.t. c6: (x[4] == -1) ==> (x[2] <= -1);
