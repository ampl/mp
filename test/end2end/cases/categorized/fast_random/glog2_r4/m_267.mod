var x{1..5} integer >= 0 <= 3;
maximize o: (-2000)*x[1] + (1000)*x[2] + (2000)*x[3] + (34)*x[4] + (3000)*x[5];
s.t. c1: abs(x[2] - x[4]) == 2;
s.t. c2: (if x[1] <= 1 then x[2] else x[5]) != 1;
