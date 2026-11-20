var x{1..4} integer >= -2 <= 0;

maximize o: 6*x[1] - 4*x[2] + 6*x[3] - 2*x[4];

s.t. c1: numberof -3 in (x[1], x[2], x[3], x[4]) <= 3;
s.t. c2: atleast 3 {i in 1..4} (x[i] >= -3);
s.t. c3: abs(x[4] - x[4]) <= 1;
