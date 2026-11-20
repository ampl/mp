var x{1..5} integer >= -3 <= 1;

minimize o: -2*x[1] - 2*x[2] + 3*x[3] - x[4] + 3*x[5];

s.t. c1: abs(x[5] - x[2]) == 0;
s.t. c2: abs(x[3] - x[5]) != 3;
s.t. c3: exists {i in 1..5} x[i] == 0;
s.t. c4: (x[1] > 1) and (x[3] <= 1) ==> (x[2] >= 1);
