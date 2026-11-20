var x{1..3} integer >= 0 <= 2;

minimize o: -x[1] + 4*x[2] + 0*x[3];

s.t. c1: (max(x[3], x[2]) > 2) <==> (x[1] == 0);
s.t. c2: (max(x[1], x[3]) <= 0) ==> (max(x[1], x[3]) <= 0);
