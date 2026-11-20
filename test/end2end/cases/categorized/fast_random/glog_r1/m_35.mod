var x{1..3} integer >= 0 <= 3;

minimize o: 5*x[1] + 6*x[2] + 5*x[3];

s.t. c1: (if x[2] >= 3 then x[3] else x[2]) > 0;
s.t. c2: (x[1] < 1) ==> (x[3] > 3) else (x[1] <= 0);
s.t. c3: ((if x[3] >= 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
