var x{1..3} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (-3)*x[2] + (-1)*x[3];
s.t. c1: (if x[3] > 1 then x[3] else x[1]) == 1;
s.t. c2: (if x[3] > 1 then x[3] else x[1]) == 1;
