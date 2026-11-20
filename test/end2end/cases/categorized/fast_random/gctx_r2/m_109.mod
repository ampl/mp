var x{1..5} integer >= 0 <= 3;

minimize o: x[1] - 2*x[2] - x[3] - 4*x[4] + 2*x[5];

s.t. c1: (x[1] >= 0) and ((x[4] < 2) or not (x[1] >= 0));
s.t. c2: (x[4] < 2) or (x[1] >= 0);
s.t. c3: (x[4] < 2) <==> (x[1] >= 0);
