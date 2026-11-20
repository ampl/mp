var x{1..4} integer >= -2 <= 1;
maximize o: (1)*x[1] + (-5)*x[2] + (-1)*x[3] + (-5)*x[4];
s.t. c1: alldiff(x[4], x[3], x[2], x[1]);
s.t. c2: (x[2] != 1) ==> (x[2] == 0) else (x[4] == -2);
s.t. c3: abs(x[4] - x[3]) == 1;
s.t. c4: atleast 2 {i in 1..4} (x[i] >= -2);
