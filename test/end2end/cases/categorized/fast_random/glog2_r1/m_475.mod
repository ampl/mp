var x{1..4} integer >= 0 <= 3;

minimize o: 3*x[1] + 2*x[2] - 4*x[3] + 0*x[4];

s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4]) != 1;
s.t. c2: exists {i in 1..4} x[i] == 1;
s.t. c3: (if x[2] > 2 then x[2] else x[2]) != 1;
s.t. c4: (if x[2] <= 2 then x[4] else x[1]) > 2;
s.t. c5: numberof 0 in (x[1], x[2], x[3], x[4]) < 0;
s.t. c6: (x[4] == 0) and (x[4] <= 0) ==> (x[4] >= 0);
s.t. c7: (x[2] == 3) ==> (x[3] != 3);
