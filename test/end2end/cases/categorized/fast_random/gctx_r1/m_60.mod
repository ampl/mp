var x{1..3} integer >= 0 <= 2;
minimize o: (-2)*x[1] + (-3)*x[2] + (6)*x[3];
s.t. c1: (max(x[2], x[2]) != 0) or (min(x[2], x[3]) < 2);
s.t. c2: not ((min(x[2], x[3]) < 2) and (min(x[2], x[3]) < 2));
s.t. c3: (max(x[2], x[2]) != 0) ==> ((max(x[2], x[2]) != 0) or (max(x[2], x[2]) != 0));
s.t. c4: not ((max(x[2], x[2]) != 0) and (max(x[2], x[2]) != 0));
