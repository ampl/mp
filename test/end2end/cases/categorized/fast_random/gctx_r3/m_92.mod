var x{1..4} integer >= -1 <= 2;

maximize o: -3*x[1] - 3*x[2] + 4*x[3] - 4*x[4];

s.t. c1: not (x[1] < -1);
s.t. c2: (if abs(x[4] - x[3]) != 3 then x[1] else x[4]) == 2;
s.t. c3: (x[1] < -1) <==> (abs(x[4] - x[3]) != 3);
