var x{1..4} integer >= -2 <= 2;

minimize o: -5*x[1] - 2*x[2] + 6*x[3] + 0*x[4];

s.t. c1: (if x[1] != -1 then x[1] else x[1]) == -1;
s.t. c2: (x[2] != 2) ==> (x[2] >= -2) else (x[3] <= -1);
