var x{1..3} integer >= 0 <= 3;

minimize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: ((if x[3] == 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c2: (x[1] < 3) or (x[1] == 3) or (x[2] < 1);
s.t. c3: (x[2] < 1) and (x[3] >= 1) ==> (x[3] >= 1);
