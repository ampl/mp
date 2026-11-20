var x{1..5} integer >= 0 <= 3;

minimize o: -2*x[1] - 5*x[2] - 4*x[3] - 5*x[4] - 4*x[5];

s.t. c1: (max(x[2], x[3]) >= 0) or (x[4] >= 0);
s.t. c2: (x[4] >= 0) <==> (max(x[2], x[3]) >= 0);
s.t. c3: (x[4] >= 0) and ((max(x[2], x[3]) >= 0) or not (x[4] >= 0));
s.t. c4: (if x[4] >= 0 then x[1] else x[1]) <= 2;
