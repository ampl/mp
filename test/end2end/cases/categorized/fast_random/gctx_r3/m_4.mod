var x{1..5} integer >= 0 <= 3;
minimize o: (1)*x[1] + (6)*x[2] + (-4)*x[3] + (5)*x[4] + (1)*x[5];
s.t. c1: (if x[5] + x[2] >= 2 then x[4] else x[3]) < 0;
s.t. c2: (max(x[1], x[5]) != 2) ==> ((max(x[3], x[1]) < 0) or (max(x[3], x[1]) < 0));
s.t. c3: (max(x[3], x[1]) < 0) or (max(x[1], x[5]) != 2);
