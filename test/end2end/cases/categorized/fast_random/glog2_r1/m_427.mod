var x{1..4} integer >= 0 <= 4;
minimize o: (-333)*x[1] + (-3000)*x[2] + (17)*x[3] + (2000)*x[4];
s.t. c1: max(x[4], x[4]) - min(x[2], x[4]) != 4;
