var x{1..3} integer >= 0 <= 3;

minimize o: -x[1] - 5*x[2] - 2*x[3];

s.t. c1: (min(x[2], x[1]) >= 3) <==> (x[2] <= 0);
s.t. c2: (x[2] <= 0) and ((min(x[2], x[1]) >= 3) or not (x[2] <= 0));
s.t. c3: (min(x[2], x[1]) >= 3) ==> ((min(x[2], x[1]) >= 3) or (x[2] <= 0));
