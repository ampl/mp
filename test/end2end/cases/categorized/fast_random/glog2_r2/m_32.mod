var x{1..5} integer >= -2 <= 0;
maximize o: (-3)*x[1] + (-5)*x[2] + (-1)*x[3] + (3)*x[4] + (-5)*x[5];
s.t. c1: abs(x[2] - x[3]) > 0;
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) != 2;
s.t. c3: (if x[4] == 0 then x[4] else x[4]) < -1;
s.t. c4: exists {i in 1..5} x[i] == 0;
