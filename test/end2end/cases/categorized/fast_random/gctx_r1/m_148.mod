var x{1..5} integer >= -1 <= 2;
maximize o: (0)*x[1] + (-2)*x[2] + (-2)*x[3] + (5)*x[4] + (-2)*x[5];
s.t. c1: not ((x[5] >= x[3]) and (x[5] >= x[3]));
s.t. c2: (if x[3] + x[5] < -1 then x[4] else x[4]) < 0;
s.t. c3: (x[5] >= x[3]) or (x[3] + x[5] < -1);
