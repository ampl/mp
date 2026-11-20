var x{1..5} integer >= 0 <= 4;
minimize o: (3000)*x[1] + (-999)*x[2] + (2000)*x[3] + (-666)*x[4] + (-333)*x[5];
s.t. c1: (if x[2] <= 0 then x[3] else x[4]) != 1;
s.t. c2: exists {i in 1..5} x[i] == 2;
