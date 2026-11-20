var x{1..3} integer >= 0 <= 3;

minimize o: -4*x[1] + 5*x[2] + 5*x[3];

s.t. c1: not ((x[3] != 1) and (x[1] > 3));
s.t. c2: (x[1] == 1) <==> (x[1] + x[1] < 1);
