var x{1..3} integer >= 0 <= 2;

maximize o: -5*x[1] + 6*x[2] - 5*x[3];

s.t. c1: (abs(x[1] - x[2]) > 2) <==> (abs(x[1] - x[2]) > 2);
s.t. c2: (max(x[1], x[3]) < 2) or (max(x[1], x[3]) < 2);
