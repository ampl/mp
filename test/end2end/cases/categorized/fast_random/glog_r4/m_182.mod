var x{1..5} integer >= -2 <= 1;
maximize o: (-5)*x[1] + (4)*x[2] + (5)*x[3] + (-3)*x[4] + (2)*x[5];
s.t. c1: (if x[3] <= -1 then x[2] else x[5]) >= -2;
s.t. c2: exists {i in 1..5} x[i] == -2;
