var x{1..4} integer >= 0 <= 3;
maximize o: (2)*x[1] + (3)*x[2] + (3)*x[3] + (-1)*x[4];
s.t. c1: exists {i in 1..4} x[i] == 2;
s.t. c2: ((x[2] < 3) ==> (x[2] <= 2)) and ((x[2] <= 2) ==> (x[2] >= 3));
s.t. c3: alldiff(x[1], x[3]);
