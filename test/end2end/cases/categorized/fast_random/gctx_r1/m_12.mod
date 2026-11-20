var x{1..3} integer >= -1 <= 2;

minimize o: -4*x[1] + 2*x[2] + 0*x[3];

s.t. c1: (x[2] < 2) ==> (x[2] < 2);
s.t. c2: (min(x[3], x[1]) > 1) <==> (abs(x[2] - x[2]) > 0);
