var x{1..5} integer >= 0 <= 4;
maximize o: (-1)*x[1] + (6)*x[2] + (1)*x[3] + (-1)*x[4] + (-1)*x[5];
s.t. c1: (x[1] == 0) ==> (x[5] < -1);
s.t. c2: alldiff(x[3], x[1]);
s.t. c3: (if x[5] >= 4 then x[4] else x[2]) != 4;
