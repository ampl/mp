var x{1..3} integer >= -2 <= 0;

minimize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: (x[3] != 0) and (x[2] == -2) ==> (x[1] >= 0);
