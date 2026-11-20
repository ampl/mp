var x{1..3} integer >= 0 <= 3;
maximize o: (-999)*x[1] + (1000)*x[2] + (51)*x[3];
s.t. c1: max(x[1], x[3]) - min(x[3], x[1]) <= 0;
s.t. c2: alldiff(x[1], x[3]);
s.t. c3: not ((x[3] >= 0) and (x[3] == 3));
