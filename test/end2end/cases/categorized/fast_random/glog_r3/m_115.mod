var x{1..6} integer >= -2 <= 0;
maximize o: (-4)*x[1] + (3)*x[2] + (4)*x[3] + (-1)*x[4] + (-1)*x[5] + (2)*x[6];
s.t. c1: exists {i in 1..6} x[i] == -2;
s.t. c2: atleast 3 {i in 1..6} (x[i] >= -2);
s.t. c3: not ((x[1] < -1) and (x[5] > -1));
s.t. c4: (if x[4] == 0 then x[2] else x[1]) < -2;
