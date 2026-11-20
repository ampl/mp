var x{1..3} integer >= -2 <= 1;

minimize o: 4*x[1] + 2*x[2] + 4*x[3];

s.t. c1: exists {i in 1..3} x[i] == -1;
s.t. c2: (if x[2] > -2 then x[2] else x[3]) == -2;
s.t. c3: (if x[3] > 0 then x[1] else x[3]) >= 0;
s.t. c4: exists {i in 1..3} x[i] == -1;
s.t. c5: (if x[3] > 0 then x[1] else x[3]) >= 0;
