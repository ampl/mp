var x{1..3} integer >= 0 <= 3;
maximize o: (2)*x[1] + (2)*x[2] + (1)*x[3];
s.t. c1: abs(x[3] - x[2]) > 2;
s.t. c2: exists {i in 1..3} x[i] == 0;
s.t. c3: alldiff(x[3], x[2]);
s.t. c4: (x[2] == 3) <==> (x[1] + x[1] < 3);
s.t. c5: (x[2] > 2) or (x[2] >= 2) or (x[3] != 2);
s.t. c6: not ((x[1] == 1) and (x[1] == 0));
s.t. c7: (if x[3] >= 0 then x[3] else x[3]) >= 2;
