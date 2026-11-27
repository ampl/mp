var x{1..3} integer >= 0 <= 3;

maximize o: -x[1] + 4*x[2] - 3*x[3];

s.t. c1: (abs(x[1] - x[1]) != 2) <==> (abs(x[1] - x[1]) != 2);
s.t. c2: (min(x[3], x[2]) < 2) <==> (abs(x[1] - x[1]) != 2);
