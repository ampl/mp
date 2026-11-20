var x{1..3} integer >= 0 <= 3;
minimize o: (5)*x[1] + (4)*x[2] + (5)*x[3];
s.t. c1: (x[2] != 0) <==> (x[1] <= 3);
s.t. c2: (if x[2] != 0 then x[2] else 0) + (if x[2] != 0 then x[2] else 0) >= 0;
