var x{1..4} integer >= 0 <= 4;

minimize o: -5*x[1] - 2*x[2] + 5*x[3] + 3*x[4];

s.t. c1: ((if x[4] == 4 then 1 else 0) + (if x[1] >= 4 then 1 else 0)) >= 1;
s.t. c2: (x[3] >= 2) ==> (x[4] != 3) else (x[3] == 0);
s.t. c3: ((x[1] <= 0) ==> (x[3] <= 3)) and ((x[3] <= 3) ==> (x[2] >= 0));
