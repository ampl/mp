var x{1..5} integer >= -2 <= 2;

maximize o: 2*x[1] - 5*x[2] - 4*x[3] - x[4] - 5*x[5];

s.t. c1: (x[3] != 2) ==> (x[1] < -1);
s.t. c2: max(x[4], x[4]) - min(x[1], x[3]) <= 2;
s.t. c3: max(x[5], x[4]) - min(x[2], x[2]) > 1;
s.t. c4: alldiff(x[3], x[5]);
s.t. c5: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) >= 1;
s.t. c6: (x[5] == -2) <==> (x[3] + x[1] == -2);
