var x{1..3} integer >= 0 <= 3;

minimize o: -5*x[1] - x[2] + 4*x[3];

s.t. c1: (if x[2] >= 0 then x[3] else x[1]) <= 0;
s.t. c2: exists {i in 1..3} x[i] == 1;
s.t. c3: x[3] == 2;
s.t. c4: x[1] == 2;
