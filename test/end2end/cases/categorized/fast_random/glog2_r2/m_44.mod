var x{1..4} integer >= -2 <= 3;

minimize o: -4*x[1] - 2*x[2] + 5*x[3] + 3*x[4];

s.t. c1: exists {i in 1..4} x[i] == 2;
s.t. c2: ((x[3] != -2) ==> (x[2] >= 2)) and ((x[2] >= 2) ==> (x[2] >= -2));
