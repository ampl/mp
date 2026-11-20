var x{1..3} integer >= -2 <= 2;
minimize o: (-5)*x[1] + (1)*x[2] + (3)*x[3];
s.t. c1: not ((x[3] != -1) and (x[1] < 2));
s.t. c2: (if x[2] <= 0 then x[1] else x[2]) >= -1;
s.t. c3: (x[3] <= 2) ==> (x[3] == -1) else (x[1] <= 2);
