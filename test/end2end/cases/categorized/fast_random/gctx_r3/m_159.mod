var x{1..5} integer >= 0 <= 2;
minimize o: (1)*x[1] + (0)*x[2] + (4)*x[3] + (-5)*x[4] + (1)*x[5];
s.t. c1: (min(x[4], x[4]) < 0) <==> (min(x[4], x[4]) < 0);
s.t. c2: (min(x[3], x[1]) != 2) and ((min(x[3], x[1]) != 2) or not (min(x[4], x[4]) < 0));
s.t. c3: (if min(x[3], x[1]) != 2 then x[1] else x[3]) == 1;
s.t. c4: (min(x[4], x[4]) < 0) ==> (min(x[4], x[4]) < 0);
