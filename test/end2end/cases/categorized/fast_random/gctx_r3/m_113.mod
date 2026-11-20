var x{1..3} integer >= 0 <= 2;

maximize o: 6*x[1] + x[2] - 5*x[3];

s.t. c1: (x[1] > x[3]) <==> (abs(x[2] - x[1]) < 2);
s.t. c2: (abs(x[2] - x[1]) < 2) ==> (abs(x[2] - x[1]) < 2);
