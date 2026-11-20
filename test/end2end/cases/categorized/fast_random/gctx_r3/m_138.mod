var x{1..3} integer >= -1 <= 1;

minimize o: 3*x[1] - 4*x[2] + 4*x[3];

s.t. c1: (x[1] == x[2]) and ((x[1] + x[2] < -1) or not (x[1] + x[2] < -1));
s.t. c2: (if x[1] == x[2] then x[3] else x[3]) == 0;
s.t. c3: (x[3] + x[1] > 1) <==> (x[3] + x[1] > 1);
s.t. c4: not (x[3] + x[1] > 1);
