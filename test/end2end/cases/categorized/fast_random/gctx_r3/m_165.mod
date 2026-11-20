var x{1..5} integer >= -1 <= 1;

maximize o: 5*x[1] + 2*x[2] - 4*x[3] + x[4] + 6*x[5];

s.t. c1: not (max(x[1], x[3]) == 0);
s.t. c2: (x[3] > x[1]) or (x[3] > x[1]);
s.t. c3: (x[3] > x[1]) and ((x[3] > x[1]) or not (x[3] > x[1]));
