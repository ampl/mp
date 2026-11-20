var x{1..3} integer >= 0 <= 2;
minimize o: (-4)*x[1] + (-1)*x[2] + (5)*x[3];
s.t. c1: alldiff(x[1], x[2]);
s.t. c2: numberof -1 in (x[1], x[2], x[3]) <= 1;
