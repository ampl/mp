var x{1..3} integer >= 0 <= 4;

minimize o: 6*x[1] + 5*x[2] - x[3];

s.t. c1: (x[1] <= 4) and (x[3] > 5) ==> (x[3] >= 4);
