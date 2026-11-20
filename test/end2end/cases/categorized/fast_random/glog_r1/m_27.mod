var x{1..5} integer >= -2 <= 2;
minimize o: (-1)*x[1] + (0)*x[2] + (4)*x[3] + (-1)*x[4] + (3)*x[5];
s.t. c1: alldiff(x[4], x[2]);
