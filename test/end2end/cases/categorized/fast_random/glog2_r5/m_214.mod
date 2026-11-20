var x{1..3} integer >= -3 <= 1;

minimize o: 4*x[1] - 3*x[2] + 3*x[3];

s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: exists {i in 1..3} x[i] == 1;
s.t. c3: (if x[2] < 0 then x[3] else x[3]) < -2;
s.t. c4: exists {i in 1..3} x[i] == 1;
