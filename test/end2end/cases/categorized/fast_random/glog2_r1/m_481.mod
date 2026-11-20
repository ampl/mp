var x{1..3} integer >= 0 <= 4;

minimize o: x[1] + 3*x[2] - 2*x[3];

s.t. c1: (x[1] == -1) ==> (x[2] < 0);
s.t. c2: abs(x[1] - x[2]) < 0;
s.t. c3: ((if x[3] > 4 then 1 else 0) + (if x[3] >= 4 then 1 else 0)) >= 1;
s.t. c4: atmost 2 {i in 1..3} (x[i] >= -1);
