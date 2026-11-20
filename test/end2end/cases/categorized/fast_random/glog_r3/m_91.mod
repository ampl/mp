var x{1..3} integer >= 0 <= 2;

minimize o: 5*x[1] - 3*x[2] + 5*x[3];

s.t. c1: (x[1] < 0) ==> (x[3] <= 2);
s.t. c2: (x[1] == 1) ==> (x[3] > 0) else (x[3] <= 0);
s.t. c3: (if x[1] <= 2 then x[3] else x[3]) <= 1;
