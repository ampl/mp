var x{1..4} integer >= 0 <= 4;

minimize o: 6*x[1] - x[2] + 2*x[3] + 3*x[4];

s.t. c1: exists {i in 1..4} x[i] == 0;
s.t. c2: ((if x[2] < 4 then 1 else 0) + (if x[1] >= 4 then 1 else 0)) >= 1;
s.t. c3: not ((x[4] > 4) and (x[1] == 0));
s.t. c4: (x[2] > 0) ==> (x[3] > 1);
