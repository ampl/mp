var x{1..4} integer >= 0 <= 4;

maximize o: 2*x[1] - 4*x[2] + 3*x[3] + 3*x[4];

s.t. c1: alldiff(x[2], x[4], x[3], x[1]);
s.t. c2: abs(x[3] - x[3]) != 4;
s.t. c3: numberof 2 in (x[1], x[2], x[3], x[4]) != 3;
s.t. c4: numberof 2 in (x[1], x[2], x[3], x[4]) < 0;
s.t. c5: abs(x[3] - x[4]) <= 2;
s.t. c6: (x[2] > 4) <==> (x[2] + x[4] == 4);
