var x{1..3} integer >= 0 <= 4;
maximize o: (2)*x[1] + (2)*x[2] + (6)*x[3];
s.t. c1: exists {i in 1..3} x[i] == 4;
s.t. c2: alldiff(x[1], x[2]);
s.t. c3: (x[3] < 1) <==> (x[2] + x[1] < 2);
s.t. c4: abs(x[2] - x[3]) < 1;
s.t. c5: exists {i in 1..3} x[i] == 4;
s.t. c6: abs(x[2] - x[3]) < 1;
