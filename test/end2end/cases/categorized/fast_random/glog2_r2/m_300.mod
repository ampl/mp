var x{1..3} integer >= 0 <= 2;
maximize o: (5)*x[1] + (4)*x[2] + (3)*x[3];
s.t. c1: (x[2] == 2) and (x[2] < 0) ==> (x[1] >= 2);
s.t. c2: atleast 3 {i in 1..3} (x[i] >= 1);
s.t. c3: alldiff(x[3], x[1]);
s.t. c4: (x[3] > 2) <==> (x[2] + x[2] != 2);
