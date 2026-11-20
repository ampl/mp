var x{1..3} integer >= 0 <= 3;

minimize o: -2*x[1] + 4*x[2] + x[3];

s.t. c1: (x[2] < 0) <==> (x[2] + x[2] >= 1);
s.t. c2: abs(x[1] - x[3]) > 2;
