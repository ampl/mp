var x{1..3} integer >= -1 <= 2;
minimize o: (2)*x[1] + (-1)*x[2] + (2)*x[3];
s.t. c1: (x[3] == 1) ==> ((min(x[2], x[1]) != 1) or (min(x[2], x[1]) != 1));
s.t. c2: (x[3] == 1) and ((x[2] == -1) or not (x[2] == -1));
