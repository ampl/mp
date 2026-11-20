var x{1..5} integer >= 0 <= 2;

minimize o: 0*x[1] - x[2] + 6*x[3] + 4*x[4] - 4*x[5];

s.t. c1: (x[5] <= 2) <==> (x[5] <= 2);
s.t. c2: (if x[5] <= 2 then x[5] else x[5]) == 1;
s.t. c3: (min(x[1], x[2]) > 0) ==> (x[3] != 1) else (x[3] != 1);
