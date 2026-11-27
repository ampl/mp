var x{1..3} integer >= 0 <= 4;

minimize o: -2*x[1] - x[2] - 3*x[3];

s.t. c1: (x[2] != 0) ==> (x[1] == 2) else (x[1] != 2);
s.t. c2: ((if x[3] >= 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
