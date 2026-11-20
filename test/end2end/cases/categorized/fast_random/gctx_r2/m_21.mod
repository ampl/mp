var x{1..3} integer >= -1 <= 1;
minimize o: (-1)*x[1] + (-4)*x[2] + (5)*x[3] + (1)*(if min(x[3], x[1]) < 1 then 1 else 0);
s.t. c1: (min(x[3], x[1]) < 1) and ((x[2] < -1) or not (x[2] < -1));
s.t. c2: (x[2] < -1) ==> (min(x[1], x[3]) == 1) else (x[2] < -1);
s.t. c3: (min(x[3], x[1]) < 1) <==> (min(x[1], x[3]) == 1);
s.t. c4: (min(x[3], x[1]) < 1) and ((x[2] < -1) or not (min(x[3], x[1]) < 1));
