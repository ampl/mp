var x{1..3} integer >= -2 <= 2;
minimize o: (-2)*x[1] + (2)*x[2] + (-3)*x[3];
s.t. c1: max(x[2], x[1]) - min(x[3], x[3]) <= -1;
s.t. c2: (if x[1] != 2 then x[1] else x[1]) < 0;
s.t. c3: ((x[3] < 1) ==> (x[3] != -2)) and ((x[3] != -2) ==> (x[3] >= 1));
