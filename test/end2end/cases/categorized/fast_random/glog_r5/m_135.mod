var x{1..3} integer >= -2 <= 1;

minimize o: x[1] + x[2] + 6*x[3];

s.t. c1: (x[3] < 0) ==> (x[1] > -2);
s.t. c2: atmost 1 {i in 1..3} (x[i] >= 0);
s.t. c3: ((if x[2] >= -1 then 1 else 0) + (if x[2] >= -1 then 1 else 0)) >= 1;
