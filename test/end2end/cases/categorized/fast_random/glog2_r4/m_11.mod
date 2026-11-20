var x{1..5} integer >= 0 <= 2;

maximize o: 6*x[1] - x[2] + 2*x[3] + 5*x[4] - x[5];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) < 0;
s.t. c2: numberof -1 in (x[1], x[2], x[3], x[4], x[5]) == 2;
s.t. c3: abs(x[5] - x[5]) == 0;
s.t. c4: (x[2] >= -1) <==> (x[4] + x[2] == 2);
