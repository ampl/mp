var x{1..3} integer >= 0 <= 2;

minimize o: 3*x[1] + 4*x[2] + 4*x[3];

s.t. c1: (x[2] > 1) <==> (x[1] + x[2] <= 1);
s.t. c2: (if x[2] == 2 then x[2] else x[2]) >= 2;
s.t. c3: (x[1] <= 2) or (x[2] >= 0) or (x[1] == 1);
s.t. c4: x[3] == 2;
