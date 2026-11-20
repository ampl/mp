var x{1..6} integer >= 0 <= 3;

maximize o: -4*x[1] + 2*x[2] + 5*x[3] + x[4] - 3*x[5] - 5*x[6];

s.t. c1: atleast 3 {i in 1..6} (x[i] >= 3);
s.t. c2: (if x[4] >= 1 then x[3] else x[6]) > 3;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 2;
s.t. c4: (if x[4] != 1 then x[3] else x[5]) < 1;
s.t. c5: x[3] == 1;
