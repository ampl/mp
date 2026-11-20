var x{1..3} integer >= -2 <= 0;

minimize o: 5*x[1] + 0*x[2] + 4*x[3];

s.t. c1: (if x[3] == 0 then x[3] else x[1]) > 0;
s.t. c2: abs(x[2] - x[1]) <= 2;
s.t. c3: exists {i in 1..3} x[i] == -2;
s.t. c4: (if x[3] != 0 then x[3] else x[1]) == -1;
s.t. c5: forall {i in 1..2} x[i] != x[i+1];
s.t. c6: (x[1] > 0) ==> (x[3] < -2) else (x[2] != 0);
