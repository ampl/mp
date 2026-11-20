var x{1..3} integer >= -2 <= 0;

minimize o: -4*x[1] + 0*x[2] - 2*x[3];

s.t. c1: (x[2] <= -2) <==> (x[2] + x[1] < 0);
s.t. c2: not ((x[1] <= -1) and (x[1] > -2));
