var x{1..3} integer >= 0 <= 2;

minimize o: 6*x[1] - x[2] - x[3];

s.t. c1: not ((abs(x[2] - x[3]) != 0) and (x[1] >= x[2]));
s.t. c2: (x[1] >= x[2]) ==> (abs(x[1] - x[3]) > 0) else (x[1] >= x[2]);
