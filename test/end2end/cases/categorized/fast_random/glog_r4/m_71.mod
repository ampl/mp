var x{1..6} integer >= 0 <= 4;
maximize o: (3)*x[1] + (2)*x[2] + (-2)*x[3] + (5)*x[4] + (6)*x[5] + (6)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 3;
s.t. c2: atleast 1 {i in 1..6} (x[i] >= 3);
s.t. c3: (if x[4] < 2 then x[2] else x[5]) > 1;
