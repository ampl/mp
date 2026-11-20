var x{1..3} integer >= 0 <= 2;

minimize o: 5*x[1] + 2*x[2] - 5*x[3];

s.t. c1: (x[1] > x[2]) ==> (max(x[3], x[3]) == 1) else (x[1] + x[1] == 1);
s.t. c2: not ((max(x[3], x[3]) == 1) and (x[1] > x[2]));
