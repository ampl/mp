var x{1..6} integer >= -2 <= 1;

maximize o: -3*x[1] + 3*x[2] - 3*x[3] - 4*x[4] - 3*x[5] + 0*x[6];

s.t. c1: max(x[1], x[4]) - min(x[3], x[5]) == 1;
s.t. c2: ((if x[1] < -2 then 1 else 0) + (if x[1] >= -2 then 1 else 0)) >= 1;
s.t. c3: exists {i in 1..6} x[i] == -2;
s.t. c4: exists {i in 1..6} x[i] == 1;
s.t. c5: (x[5] != 0) ==> (x[2] <= 1);
