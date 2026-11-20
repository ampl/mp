var x{1..3} integer >= -1 <= 1;

minimize o: 0*x[1] + 3*x[2] + 3*x[3];

s.t. c1: (max(x[2], x[2]) == -1) ==> (x[1] != x[2]);
s.t. c2: (min(x[1], x[2]) == 1) ==> ((x[1] != x[2]) or (max(x[2], x[2]) == -1));
