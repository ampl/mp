var x{1..5} integer >= 0 <= 1;

minimize o: -x[1] + 6*x[2] + 2*x[3] + 0*x[4] - 5*x[5];

s.t. c1: exists {i in 1..5} x[i] == 0;
s.t. c2: atmost 3 {i in 1..5} (x[i] >= 0);
s.t. c3: ((if x[5] < 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c4: (x[3] > 0) or (x[1] >= 1) or (x[5] <= 1);
