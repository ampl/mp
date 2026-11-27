var x{1..5} integer >= -2 <= 2;

maximize o: 4*x[1] + 0*x[2] - x[3] - 3*x[4] + x[5];

s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) == 0;
s.t. c4: (if x[5] != 2 then x[4] else x[1]) != 2;
