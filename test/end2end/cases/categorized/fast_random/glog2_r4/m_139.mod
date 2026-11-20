var x{1..4} integer >= 0 <= 4;
minimize o: (-5)*x[1] + (2)*x[2] + (-5)*x[3] + (6)*x[4];
s.t. c1: alldiff(x[3], x[1]);
s.t. c2: (x[2] == 2) and (x[3] == 4) ==> (x[4] >= 2);
s.t. c3: alldiff(x[3], x[1]);
s.t. c4: (x[2] == 2) and (x[3] == 4) ==> (x[4] >= 2);
