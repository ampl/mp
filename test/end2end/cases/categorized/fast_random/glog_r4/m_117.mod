var x{1..3} integer >= 0 <= 4;
minimize o: (2)*x[1] + (3)*x[2] + (2)*x[3];
s.t. c1: ((x[2] != 3) ==> (x[1] == 4)) and ((x[1] == 4) ==> (x[3] >= 3));
s.t. c2: alldiff(x[3], x[1]);
