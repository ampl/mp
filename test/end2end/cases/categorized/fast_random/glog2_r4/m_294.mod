var x{1..5} integer >= -3 <= 1;
minimize o: (1)*x[1] + (4)*x[2] + (3)*x[3] + (4)*x[4] + (-3)*x[5];
s.t. c1: (x[1] < -2) ==> (x[3] > 1) else (x[5] == -2);
s.t. c2: alldiff(x[3], x[1], x[5], x[4]);
s.t. c3: alldiff(x[5], x[1]);
