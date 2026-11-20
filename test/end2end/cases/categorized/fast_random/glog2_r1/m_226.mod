var x{1..4} integer >= 0 <= 2;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4];

s.t. c1: (x[2] >= 0) ==> (x[3] == 1);
s.t. c2: ((if x[2] < 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
