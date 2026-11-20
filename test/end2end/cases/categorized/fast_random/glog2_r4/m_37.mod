var x{1..3} integer >= 0 <= 2;

minimize o: 34*x[1] - 666*x[2] + 34*x[3];

s.t. c1: (x[2] > 1) and (x[2] != 0) ==> (x[2] >= 1);
