var x{1..3} integer >= 0 <= 4;
minimize o: (-2000)*x[1] + (-999)*x[2] + (34)*x[3];
s.t. c1: max(x[1], x[2]) - min(x[2], x[2]) > 0;
s.t. c2: alldiff(x[1], x[3]);
