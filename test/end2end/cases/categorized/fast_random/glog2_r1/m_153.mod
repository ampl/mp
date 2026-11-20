var x{1..3} integer >= -2 <= 1;

minimize o: -4*x[1] - 3*x[2] + 6*x[3];

s.t. c1: (x[1] < 1) and (x[2] > -2) ==> (x[1] >= 1);
s.t. c2: x[2] == 1;
