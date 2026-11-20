var x{1..6} integer >= 0 <= 3;
maximize o: (2)*x[1] + (-2)*x[2] + (-4)*x[3] + (-1)*x[4] + (-5)*x[5] + (2)*x[6];
s.t. c1: not ((x[2] > 1) and (x[5] <= 3));
s.t. c2: atleast 2 {i in 1..6} (x[i] >= 1);
s.t. c3: alldiff(x[5], x[1], x[4], x[6]);
s.t. c4: (x[2] < 0) and (x[4] <= 0) ==> (x[1] >= 0);
s.t. c5: x[1] == 2;
s.t. c6: x[4] == 3;
