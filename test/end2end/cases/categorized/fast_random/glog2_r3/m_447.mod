var x{1..4} integer >= 0 <= 2;

maximize o: 3000*x[1] - 666*x[2] + 1000*x[3] - 999*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: not ((x[4] >= 1) and (x[1] <= 0));
s.t. c3: (if x[1] < 1 then x[1] else x[4]) != 0;
