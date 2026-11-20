var x{1..6} integer >= 0 <= 4;
minimize o: (-1000)*x[1] + (-1000)*x[2] + (1000)*x[3] + (17)*x[4] + (1000)*x[5] + (-999)*x[6];
s.t. c1: abs(x[1] - x[5]) <= 2;
s.t. c2: alldiff(x[2], x[1]);
s.t. c3: (x[5] < 2) <==> (x[3] + x[3] <= 2);
s.t. c4: ((if x[6] == 4 then 1 else 0) + (if x[3] >= 4 then 1 else 0)) >= 1;
