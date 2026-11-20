var x{1..4} integer >= 0 <= 4;

minimize o: -2*x[1] + x[2] + 4*x[3] + 4*x[4];

s.t. c1: ((if x[1] <= 4 then 1 else 0) + (if x[3] >= 4 then 1 else 0)) >= 1;
s.t. c2: (x[1] >= 4) ==> (x[1] == 3);
