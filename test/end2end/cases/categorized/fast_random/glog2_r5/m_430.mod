var x{1..5} integer >= 0 <= 4;

minimize o: x[1] + 3*x[2] + 5*x[3] - x[4] - 5*x[5];

s.t. c1: ((if x[2] >= 1 then 1 else 0) + (if x[4] >= 1 then 1 else 0)) >= 1;
s.t. c2: ((x[2] != 4) ==> (x[4] >= 4)) and ((x[4] >= 4) ==> (x[1] >= 4));
s.t. c3: (x[1] == 0) ==> (x[1] >= 3) else (x[1] > 0);
