var x{1..3} integer >= 0 <= 4;
maximize o: (-2000)*x[1] + (-1000)*x[2] + (2000)*x[3];
s.t. c1: not ((x[2] != 1) and (x[2] <= 1));
s.t. c2: (x[1] == 3) or (x[1] != 0) or (x[1] <= 4);
s.t. c3: alldiff(x[3], x[1]);
s.t. c4: (x[2] <= 1) ==> (x[2] > 2);
