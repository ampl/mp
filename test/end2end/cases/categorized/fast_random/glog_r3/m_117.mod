var x{1..5} integer >= -2 <= 1;
maximize o: (2)*x[1] + (6)*x[2] + (1)*x[3] + (2)*x[4] + (0)*x[5];
s.t. c1: abs(x[4] - x[5]) <= 3;
s.t. c2: alldiff(x[5], x[1]);
s.t. c3: (x[3] < -2) ==> (x[3] <= 1);
