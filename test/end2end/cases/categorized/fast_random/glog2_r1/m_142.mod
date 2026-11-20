var x{1..6} integer >= 0 <= 4;

minimize o: -2*x[1] + 2*x[2] + x[3] + 3*x[4] - 5*x[5] + 5*x[6];

s.t. c1: (x[1] == 3) and (x[2] <= 1) ==> (x[6] >= 3);
s.t. c2: exists {i in 1..6} x[i] == 0;
s.t. c3: ((if x[1] >= 4 then 1 else 0) + (if x[1] >= 4 then 1 else 0)) >= 1;
s.t. c4: (x[2] != 1) ==> (x[1] == 3);
