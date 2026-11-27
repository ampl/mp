var x{1..6} integer >= 0 <= 3;

minimize o: -2*x[1] + 5*x[2] + 6*x[3] - x[4] - 3*x[5] - 4*x[6];

s.t. c1: (x[6] > 0) ==> (x[5] <= 3);
s.t. c2: not ((x[1] > 2) and (x[5] > 1));
s.t. c3: (x[4] < 3) ==> (x[5] < 1);
s.t. c4: atleast 1 {i in 1..6} (x[i] >= 0);
s.t. c5: (if x[1] < 1 then x[6] else x[4]) > 1;
s.t. c6: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 1;
s.t. c7: forall {i in 1..5} x[i] != x[i+1];
