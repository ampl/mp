var x{1..5} integer >= -2 <= 3;

minimize o: 5*x[1] + 4*x[2] - 3*x[3] - x[4] - 2*x[5];

s.t. c1: (x[3] < 0) <==> (x[3] + x[5] >= 2);
s.t. c2: (if x[4] < 3 then x[2] else x[4]) > -2;
s.t. c3: not ((x[2] < -2) and (x[1] < 2));
