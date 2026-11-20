var x{1..6} integer >= 0 <= 3;

minimize o: -2*x[1] + 6*x[2] + 2*x[3] + 5*x[4] + 6*x[5] - 3*x[6];

s.t. c1: atmost 3 {i in 1..6} (x[i] >= 0);
s.t. c2: ((if x[1] <= 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) >= 1;
s.t. c4: alldiff(x[2], x[6], x[5], x[3]);
s.t. c5: atmost 3 {i in 1..6} (x[i] >= 0);
s.t. c6: alldiff(x[2], x[6], x[5], x[3]);
