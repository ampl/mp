var x{1..6} integer >= 0 <= 4;
maximize o: (-2)*x[1] + (3)*x[2] + (5)*x[3] + (0)*x[4] + (-4)*x[5] + (3)*x[6];
s.t. c1: atleast 3 {i in 1..6} (x[i] >= 1);
s.t. c2: (if x[5] <= 1 then x[6] else x[1]) != 0;
s.t. c3: abs(x[4] - x[6]) > 1;
