var x{1..4} integer >= 0 <= 4;

minimize o: x[1] - 5*x[2] + x[3] - 4*x[4];

s.t. c1: (x[3] >= 0) ==> (x[1] == 1) else (x[2] < 2);
s.t. c2: ((x[2] < 2) ==> (x[4] < 3)) and ((x[4] < 3) ==> (x[4] >= 2));
s.t. c3: exists {i in 1..4} x[i] == 0;
