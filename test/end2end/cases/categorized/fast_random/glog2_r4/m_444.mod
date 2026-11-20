var x{1..5} integer >= -3 <= 2;
minimize o: (-5)*x[1] + (0)*x[2] + (4)*x[3] + (2)*x[4] + (3)*x[5];
s.t. c1: alldiff(x[3], x[1], x[5], x[4]);
s.t. c2: (if x[5] < 1 then x[3] else x[1]) == 1;
s.t. c3: (x[3] == 1) and (x[1] < -1) ==> (x[5] >= 1);
