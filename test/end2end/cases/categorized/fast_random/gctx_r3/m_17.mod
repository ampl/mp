var x{1..4} integer >= -1 <= 1;

minimize o: -2*x[1] - 4*x[2] + x[3] - 3*x[4];

s.t. c1: (x[3] >= -1) and ((x[4] < -1) or not (x[4] < -1));
s.t. c2: (x[3] >= -1) <==> (x[3] >= -1);
s.t. c3: (x[4] < -1) <==> (x[3] >= -1);
s.t. c4: (if x[3] >= -1 then x[3] else x[3]) < 1;
