var x{1..4} integer >= 0 <= 3;
minimize o: (1)*x[1] + (2)*x[2] + (1)*x[3] + (-5)*x[4];
s.t. c1: (if x[1] <= 0 then x[3] else x[4]) <= 3;
