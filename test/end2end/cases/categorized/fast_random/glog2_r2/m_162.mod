var x{1..3} integer >= 0 <= 4;

minimize o: 0*x[1] - 5*x[2] + 6*x[3];

s.t. c1: (x[2] >= 4) and (x[1] > 4) ==> (x[2] >= 4);
