var x{1..4} integer >= 0 <= 3;
minimize o: (-1)*x[1] + (-1)*x[2] + (-4)*x[3] + (-4)*x[4];
s.t. c1: (min(x[1], x[2]) >= 3) ==> (x[3] == 0);
s.t. c2: (x[3] == 0) and ((x[3] == 0) or not (min(x[1], x[2]) >= 3));
