var x{1..3} integer >= 0 <= 2;
minimize o: (-1)*x[1] + (-1)*x[2] + (1)*x[3] + (-3)*(if min(x[2], x[1]) != 1 then 1 else 0);
s.t. c1: (x[1] >= x[1]) ==> (x[1] >= x[1]) else (min(x[2], x[1]) != 1);
s.t. c2: (x[1] >= x[1]) <==> (min(x[2], x[1]) != 1);
