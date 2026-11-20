var x{1..5} integer >= -2 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (-3)*x[3] + (1)*x[4] + (-3)*x[5];
s.t. c1: alldiff(x[3], x[4]);
s.t. c2: not ((x[3] > 0) and (x[3] >= 2));
s.t. c3: (x[2] == -1) ==> (x[5] <= 1) else (x[2] > 0);
s.t. c4: x[5] == 1;
s.t. c5: x[1] == -1;
