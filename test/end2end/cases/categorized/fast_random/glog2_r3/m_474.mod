var x{1..3} integer >= -2 <= 3;

minimize o: -x[1] + 2*x[2] + 6*x[3];

s.t. c1: abs(x[2] - x[1]) != 5;
s.t. c2: ((x[2] == -1) ==> (x[3] != 0)) and ((x[3] != 0) ==> (x[3] >= -1));
s.t. c3: exists {i in 1..3} x[i] == 3;
