var x{1..5} integer >= 0 <= 2;
maximize o: (-5)*x[1] + (4)*x[2] + (4)*x[3] + (2)*x[4] + (1)*x[5];
s.t. c1: (x[3] == 0) <==> (x[5] + x[5] != 0);
s.t. c2: exists {i in 1..5} x[i] == 1;
