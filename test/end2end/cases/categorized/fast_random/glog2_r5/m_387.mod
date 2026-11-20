var x{1..6} integer >= 0 <= 2;
maximize o: (-3000)*x[1] + (3000)*x[2] + (-1000)*x[3] + (3000)*x[4] + (-666)*x[5] + (51)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 1;
s.t. c2: (if x[5] == 0 then x[2] else x[2]) > 1;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) != 2;
