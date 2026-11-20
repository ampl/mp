var x{1..6} integer >= 0 <= 2;
maximize o: (-999)*x[1] + (34)*x[2] + (1000)*x[3] + (-333)*x[4] + (-333)*x[5] + (2000)*x[6];
s.t. c1: not ((x[3] != 1) and (x[4] < 2));
s.t. c2: (if x[2] <= 2 then x[2] else x[3]) == 2;
