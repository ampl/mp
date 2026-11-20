var x{1..5} integer >= -2 <= 2;
minimize o: (3)*x[1] + (2)*x[2] + (-1)*x[3] + (-1)*x[4] + (-1)*x[5];
s.t. c1: max(x[4], x[5]) - min(x[5], x[5]) > -1;
