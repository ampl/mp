var x{1..3} integer >= 0 <= 2;

minimize o: 2*x[1] - 5*x[2] - 2*x[3];

s.t. c1: (x[2] > 0) ==> (x[3] == 0);
s.t. c2: (x[3] != 2) ==> (x[1] < 3) else (x[1] < 2);
s.t. c3: ((if x[3] > 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
