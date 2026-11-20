var x{1..3} integer >= 0 <= 4;
minimize o: (-2)*x[1] + (-2)*x[2] + (1)*x[3];
s.t. c1: alldiff(x[1], x[2]);
s.t. c2: alldiff(x[1], x[2]);
