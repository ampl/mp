var x{1..4} integer >= 0 <= 2;

maximize o: -5*x[1] + 6*x[2] + 5*x[3] - 3*x[4];

s.t. c1: abs(x[2] - x[1]) <= 0;
s.t. c2: abs(x[1] - x[1]) < 1;
s.t. c3: max(x[3], x[2]) - min(x[1], x[4]) > 2;
s.t. c4: numberof 1 in (x[1], x[2], x[3], x[4]) < 3;
s.t. c5: x[2] == 1;
s.t. c6: x[3] == 1;
