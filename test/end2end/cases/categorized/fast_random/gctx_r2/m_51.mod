var x{1..4} integer >= -1 <= 1;

minimize o: 4*x[1] + 5*x[2] + x[3] - 3*x[4] + (if x[1] > 0 then 1 else 0);

s.t. c1: (x[1] < 1) and ((x[1] > 0) or not (min(x[2], x[2]) >= 0));
s.t. c2: (min(x[2], x[2]) >= 0) <==> (x[1] < 1);
s.t. c3: (x[1] > 0) <==> (min(x[2], x[2]) >= 0);
s.t. c4: (min(x[2], x[2]) >= 0) or (min(x[2], x[2]) >= 0);
