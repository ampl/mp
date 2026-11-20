var x{1..5} integer >= 0 <= 3;

minimize o: 2*x[1] - 2*x[2] - x[3] - 2*x[4] - 4*x[5];

s.t. c1: (x[1] != 2) ==> (x[4] >= 1) else (x[5] >= 1);
s.t. c2: (if x[3] == 3 then x[4] else x[1]) != 0;
s.t. c3: exists {i in 1..5} x[i] == 0;
