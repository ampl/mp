var x{1..3} integer >= 0 <= 2;

minimize o: 0*x[1] - 3*x[2] + 5*x[3];

s.t. c1: ((if x[2] < 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c2: ((x[1] >= 2) ==> (x[3] >= 1)) and ((x[3] >= 1) ==> (x[2] >= 2));
s.t. c3: ((if x[1] > 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
