var x{1..5} integer >= -2 <= 2;

maximize o: -4*x[1] - x[2] - 5*x[3] - x[4] - x[5];

s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) > 3;
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) == 1;
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
s.t. c4: ((if x[1] > -1 then 1 else 0) + (if x[1] >= -1 then 1 else 0)) >= 1;
