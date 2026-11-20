var x{1..5} integer >= -2 <= 1;
maximize o: (-5)*x[1] + (1)*x[2] + (-2)*x[3] + (3)*x[4] + (5)*x[5];
s.t. c1: not ((x[4] != 1) and (x[5] != -2));
s.t. c2: alldiff(x[5], x[4], x[1]);
s.t. c3: (x[1] > 0) or (x[1] < -1) or (x[2] != -1);
s.t. c4: abs(x[5] - x[3]) == 0;
