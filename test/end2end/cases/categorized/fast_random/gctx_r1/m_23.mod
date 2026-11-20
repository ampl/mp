var x{1..3} integer >= -1 <= 2;

maximize o: -x[1] + 4*x[2] - 5*x[3];

s.t. c1: (x[1] + x[1] >= 1) <==> (x[1] + x[1] >= 1);
s.t. c2: (if x[2] >= -1 then x[2] else x[2]) < 0;
s.t. c3: (if x[1] + x[1] >= 1 then x[2] else x[3]) <= 0;
