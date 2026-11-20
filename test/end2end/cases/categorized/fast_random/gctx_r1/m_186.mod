var x{1..4} integer >= -1 <= 2;

minimize o: 2*x[1] - x[2] - 3*x[3] + 2*x[4] + 2*(if x[2] <= 0 then 1 else 0);

s.t. c1: (x[2] <= 0) ==> (x[1] == -1) else (x[2] <= 0);
s.t. c2: (x[2] <= 0) ==> (x[1] == -1);
s.t. c3: (if x[1] == -1 then x[4] else x[3]) < 0;
