var x{1..3} integer >= -2 <= 2;

minimize o: 5*x[1] - 4*x[2] - 4*x[3];

s.t. c1: not ((x[2] >= 1) and (x[3] >= 0));
s.t. c2: (x[1] == -1) ==> (x[1] < -1);
s.t. c3: ((if x[3] >= -1 then 1 else 0) + (if x[3] >= -1 then 1 else 0)) >= 1;
