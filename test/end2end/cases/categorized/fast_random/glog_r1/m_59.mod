var x{1..5} integer >= -2 <= 1;

minimize o: -2*x[1] - 4*x[2] - 3*x[3] + 2*x[4] - x[5];

s.t. c1: not ((x[3] >= 0) and (x[5] >= 1));
s.t. c2: ((if x[4] >= -2 then 1 else 0) + (if x[2] >= -2 then 1 else 0)) >= 1;
s.t. c3: (x[5] <= -1) ==> (x[1] != 1) else (x[1] == 1);
s.t. c4: (x[5] > 0) ==> (x[1] >= 1) else (x[3] <= 0);
