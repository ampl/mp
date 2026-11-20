var x{1..3} integer >= -1 <= 2;
minimize o: (-5)*x[1] + (-1)*x[2] + (2)*x[3] + (3)*(if x[3] == 0 then 1 else 0);
s.t. c1: (x[3] == 0) ==> (x[2] != 0);
s.t. c2: not ((min(x[3], x[1]) != -1) and (x[2] != 0));
