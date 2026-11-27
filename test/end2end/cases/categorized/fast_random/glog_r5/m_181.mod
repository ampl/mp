var x{1..5} integer >= -2 <= 2;

minimize o: 3*x[1] - 2*x[2] + 5*x[3] - 4*x[4] - 2*x[5];

s.t. c1: atmost 1 {i in 1..5} (x[i] >= -1);
s.t. c2: (x[4] < 2) ==> (x[3] == -1) else (x[3] <= -1);
s.t. c3: ((if x[3] <= -1 then 1 else 0) + (if x[5] >= -1 then 1 else 0)) >= 1;
