var x{1..3} integer >= 0 <= 3;
maximize o: (-4)*x[1] + (2)*x[2] + (1)*x[3];
s.t. c1: (x[2] > 0) or (x[3] > 1) or (x[3] == 3);
s.t. c2: alldiff(x[2], x[3]);
s.t. c3: not ((x[2] <= 0) and (x[1] < 0));
s.t. c4: (x[1] != 3) and (x[1] < 0) ==> (x[2] >= 3);
s.t. c5: (x[2] > 0) or (x[3] > 1) or (x[3] == 3);
s.t. c6: (x[1] != 3) and (x[1] < 0) ==> (x[2] >= 3);
