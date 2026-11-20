var x{1..5} integer >= 0 <= 3;
minimize o: (34)*x[1] + (-3000)*x[2] + (34)*x[3] + (17)*x[4] + (-333)*x[5];
s.t. c1: (if x[3] <= 3 then x[5] else x[1]) < 0;
s.t. c2: abs(x[2] - x[1]) <= 2;
