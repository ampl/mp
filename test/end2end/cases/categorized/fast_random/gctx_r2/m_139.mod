var x{1..5} integer >= 0 <= 2;
minimize o: (5)*x[1] + (2)*x[2] + (-3)*x[3] + (0)*x[4] + (-5)*x[5];
s.t. c1: (min(x[3], x[3]) >= 2) ==> ((x[2] <= 1) or (min(x[3], x[3]) >= 2));
s.t. c2: (if x[3] + x[4] == 0 then x[3] else x[2]) < 1;
s.t. c3: (x[2] <= 1) ==> (min(x[3], x[3]) >= 2) else (x[2] <= 1);
