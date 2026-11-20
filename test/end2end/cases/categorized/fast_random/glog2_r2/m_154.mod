var x{1..5} integer >= -3 <= 1;

minimize o: -2*x[1] + 3*x[2] - 5*x[3] + x[4] + 0*x[5];

s.t. c1: not ((x[4] != 0) and (x[5] < -3));
s.t. c2: numberof -2 in (x[1], x[2], x[3], x[4], x[5]) == 1;
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
s.t. c4: (x[5] != -2) and (x[3] == 1) ==> (x[3] >= -2);
