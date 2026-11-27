var x{1..5} integer >= 0 <= 2;

maximize o: 2*x[1] - 4*x[2] + 3*x[3] + 6*x[4] + 2*x[5];

s.t. c1: (abs(x[4] - x[5]) < 1) <==> (abs(x[4] - x[5]) < 1);
s.t. c2: (abs(x[4] - x[5]) < 1) or (min(x[2], x[3]) == 1);
