var x{1..5} integer >= 0 <= 4;
minimize o: (34)*x[1] + (-1000)*x[2] + (34)*x[3] + (-3000)*x[4] + (-999)*x[5];
s.t. c1: max(x[3], x[1]) - min(x[2], x[4]) != 0;
s.t. c2: exists {i in 1..5} x[i] == 2;
s.t. c3: ((if x[2] < 2 then 1 else 0) + (if x[4] >= 2 then 1 else 0)) >= 1;
s.t. c4: max(x[3], x[3]) - min(x[3], x[5]) >= 1;
