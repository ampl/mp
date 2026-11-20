var x{1..3} integer >= -2 <= 1;
minimize o: (-999)*x[1] + (-3000)*x[2] + (51)*x[3];
s.t. c1: max(x[3], x[1]) - min(x[1], x[1]) != 0;
