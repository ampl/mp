var x{1..6} integer >= 0 <= 2;

minimize o: 6*x[1] + 0*x[2] + 5*x[3] - 2*x[4] - 3*x[5] + 6*x[6];

s.t. c1: ((if x[4] >= 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
s.t. c2: max(x[5], x[6]) - min(x[1], x[4]) > 0;
s.t. c3: (x[1] < 1) ==> (x[2] >= 1) else (x[4] <= 2);
s.t. c4: (x[4] < 2) and (x[6] >= 1) ==> (x[3] >= 2);
s.t. c5: exists {i in 1..6} x[i] == 2;
s.t. c6: ((x[4] <= 0) ==> (x[2] <= 2)) and ((x[2] <= 2) ==> (x[4] >= 0));
