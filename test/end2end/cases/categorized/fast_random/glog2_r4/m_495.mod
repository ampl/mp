var x{1..4} integer >= 0 <= 2;
minimize o: (3)*x[1] + (-4)*x[2] + (-4)*x[3] + (5)*x[4];
s.t. c1: alldiff(x[3], x[2], x[1]);
s.t. c2: (x[2] >= 0) ==> (x[1] > 2);
s.t. c3: (x[2] >= 2) <==> (x[1] + x[1] != 1);
s.t. c4: exists {i in 1..4} x[i] == 2;
s.t. c5: numberof 0 in (x[1], x[2], x[3], x[4]) != 1;
