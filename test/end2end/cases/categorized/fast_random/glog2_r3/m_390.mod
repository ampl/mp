var x{1..3} integer >= -2 <= 2;

minimize o: -2*x[1] - x[2] + 0*x[3];

s.t. c1: (x[2] == 1) and (x[2] == 1) ==> (x[3] >= 1);
s.t. c2: numberof -1 in (x[1], x[2], x[3]) == 1;
s.t. c3: (x[1] == 2) <==> (x[1] + x[3] == -1);
