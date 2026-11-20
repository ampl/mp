var x{1..3} integer >= -2 <= 1;
minimize o: (-1)*x[1] + (3)*x[2] + (3)*x[3];
s.t. c1: (if x[2] == -1 then x[2] else x[3]) != -1;
s.t. c2: alldiff(x[2], x[3], x[1]);
