var x{1..4} integer >= -1 <= 1;

minimize o: 2*x[1] - x[2] - x[3] + 4*x[4];

s.t. c1: (x[1] >= 0) ==> ((abs(x[2] - x[2]) < 2) or (min(x[3], x[3]) > 1));
s.t. c2: (min(x[3], x[3]) > 1) <==> (x[1] >= 0);
