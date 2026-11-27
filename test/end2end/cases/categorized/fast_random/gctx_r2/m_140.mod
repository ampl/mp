var x{1..3} integer >= 0 <= 3;

minimize o: -x[1] + 0*x[2] - 2*x[3];

s.t. c1: not ((x[2] + x[2] == 0) and (abs(x[3] - x[3]) != 3));
s.t. c2: (abs(x[3] - x[3]) != 3) <==> (abs(x[3] - x[3]) != 3);
