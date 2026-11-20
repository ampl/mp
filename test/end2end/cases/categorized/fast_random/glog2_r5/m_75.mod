var x{1..3} integer >= -2 <= 1;
maximize o: (-1)*x[1] + (-5)*x[2] + (-4)*x[3];
s.t. c1: (x[3] > 1) <==> (x[3] + x[3] > -1);
s.t. c2: (x[2] != -2) ==> (x[1] != 1) else (x[3] < -1);
s.t. c3: alldiff(x[2], x[1], x[3]);
s.t. c4: exists {i in 1..3} x[i] == -2;
s.t. c5: numberof 0 in (x[1], x[2], x[3]) < 2;
s.t. c6: (x[3] <= 0) or (x[3] >= -1) or (x[3] != 0);
