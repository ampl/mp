var x{1..4} integer >= 0 <= 2;

minimize o: -4*x[1] + x[2] + 0*x[3] + 3*x[4];

s.t. c1: max(x[4], x[3]) - min(x[4], x[3]) == 0;
s.t. c2: ((x[1] <= -1) ==> (x[2] > 3)) and ((x[2] > 3) ==> (x[3] >= -1));
s.t. c3: max(x[4], x[2]) - min(x[1], x[2]) == -1;
s.t. c4: exists {i in 1..4} x[i] == 3;
