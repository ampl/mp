var x{1..3} integer >= 0 <= 3;
minimize o: (1)*x[1] + (-3)*x[2] + (-5)*x[3];
s.t. c1: alldiff(x[3], x[2]);
s.t. c2: (x[1] < 1) and (x[1] >= 2) ==> (x[2] >= 1);
s.t. c3: not ((x[2] < 3) and (x[3] == 1));
s.t. c4: alldiff(x[3], x[2], x[1]);
