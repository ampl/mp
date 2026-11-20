var x{1..3} integer >= 0 <= 2;
maximize o: (3000)*x[1] + (-999)*x[2] + (34)*x[3];
s.t. c1: alldiff(x[1], x[2]);
s.t. c2: (if x[2] == 2 then x[2] else x[2]) <= 1;
