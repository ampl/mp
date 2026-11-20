var x{1..6} integer >= -2 <= 2;

minimize o: 0*x[1] + 5*x[2] - 4*x[3] + 4*x[4] - 3*x[5] + 5*x[6];

s.t. c1: (x[6] == 1) and (x[3] < 1) ==> (x[6] >= 1);
s.t. c2: not ((x[5] <= 1) and (x[6] < 1));
s.t. c3: (if x[1] < 0 then x[4] else x[5]) > 0;
s.t. c4: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) != 2;
s.t. c5: (x[4] >= 2) or (x[5] < 0) or (x[5] > 1);
s.t. c6: forall {i in 1..5} x[i] != x[i+1];
s.t. c7: (x[4] <= 1) or (x[4] == 1) or (x[2] < 2);
