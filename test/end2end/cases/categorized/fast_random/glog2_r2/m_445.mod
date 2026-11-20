var x{1..6} integer >= 0 <= 4;

minimize o: -2*x[1] - 4*x[2] + 6*x[3] + 5*x[4] + 0*x[5] - 5*x[6];

s.t. c1: ((if x[4] != 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
s.t. c2: alldiff(x[2], x[6], x[5], x[4]);
s.t. c3: (x[6] <= 2) or (x[5] >= 3) or (x[4] == 3);
s.t. c4: ((if x[2] < 4 then 1 else 0) + (if x[6] >= 4 then 1 else 0)) >= 1;
s.t. c5: atmost 3 {i in 1..6} (x[i] >= 1);
