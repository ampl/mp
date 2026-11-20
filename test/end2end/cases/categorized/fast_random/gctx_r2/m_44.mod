var x{1..3} integer >= -1 <= 1;

maximize o: 0*x[1] - 2*x[2] + 3*x[3];

s.t. c1: (x[2] == 0) <==> (x[2] == 0);
s.t. c2: (if x[3] + x[2] >= 0 then x[2] else x[1]) < -1;
s.t. c3: not (x[2] == 0);
s.t. c4: not (x[3] + x[2] >= 0);
