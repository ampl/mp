var x{1..4} integer >= -2 <= 2;
maximize o: (1)*x[1] + (4)*x[2] + (1)*x[3] + (3)*x[4];
s.t. c1: exists {i in 1..4} x[i] == 0;
s.t. c2: (x[2] > 2) <==> (x[2] + x[4] != -2);
s.t. c3: not ((x[1] < 2) and (x[3] > -1));
s.t. c4: numberof 2 in (x[1], x[2], x[3], x[4]) < 2;
s.t. c5: (x[4] < 2) <==> (x[4] + x[2] != 2);
