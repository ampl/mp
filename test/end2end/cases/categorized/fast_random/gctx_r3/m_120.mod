var x{1..3} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (-3)*x[2] + (-3)*x[3];
s.t. c1: (if min(x[3], x[2]) != 2 then x[1] else 0) + (if x[3] <= 1 then x[1] else 0) >= 0;
s.t. c2: (if min(x[3], x[2]) != 2 then x[2] else 0) + (if x[3] <= 1 then x[3] else 0) < 1;
s.t. c3: not (x[3] <= 1);
s.t. c4: (min(x[3], x[2]) != 2) <==> (x[3] <= 1);
