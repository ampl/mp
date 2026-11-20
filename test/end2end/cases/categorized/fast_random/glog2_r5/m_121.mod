var x{1..6} integer >= 0 <= 2;

minimize o: x[1] - x[2] - x[3] - 2*x[4] + 4*x[5] - 3*x[6];

s.t. c1: atmost 3 {i in 1..6} (x[i] >= -1);
s.t. c2: max(x[1], x[2]) - min(x[4], x[2]) >= 0;
s.t. c3: abs(x[1] - x[6]) != 1;
s.t. c4: ((if x[4] <= -1 then 1 else 0) + (if x[4] >= -1 then 1 else 0)) >= 1;
