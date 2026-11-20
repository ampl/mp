var x{1..6} integer >= 0 <= 2;

minimize o: -4*x[1] + 6*x[2] + 3*x[3] + 5*x[4] + 4*x[5] - 4*x[6];

s.t. c1: ((if x[1] > 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
s.t. c2: alldiff(x[6], x[2]);
s.t. c3: (if x[2] <= 1 then x[4] else x[2]) <= 0;
s.t. c4: exists {i in 1..6} x[i] == 0;
s.t. c5: x[1] == 0;
s.t. c6: x[4] == 1;
