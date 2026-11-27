var x{1..4} integer >= -2 <= 1;

minimize o: 6*x[1] + 5*x[2] + 5*x[3] + 2*x[4];

s.t. c1: ((if x[4] >= -1 then 1 else 0) + (if x[4] >= -1 then 1 else 0)) >= 1;
s.t. c2: ((x[4] > 0) ==> (x[2] < -1)) and ((x[2] < -1) ==> (x[1] >= 0));
s.t. c3: (x[3] >= -1) ==> (x[4] < -1) else (x[1] >= 0);
