var x{1..4} integer >= 0 <= 4;

minimize o: -3000*x[1] + 17*x[2] - 666*x[3] + 34*x[4];

s.t. c1: (if x[2] <= 0 then x[1] else x[1]) == 2;
s.t. c2: (if x[4] >= 3 then x[3] else x[3]) != 1;
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
