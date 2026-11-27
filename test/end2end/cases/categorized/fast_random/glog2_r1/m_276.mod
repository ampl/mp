var x{1..4} integer >= 0 <= 4;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4];

s.t. c1: atleast 2 {i in 1..4} (x[i] >= 3);
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4]) < 0;
s.t. c3: (x[4] < 2) or (x[1] <= 4) or (x[2] > 0);
