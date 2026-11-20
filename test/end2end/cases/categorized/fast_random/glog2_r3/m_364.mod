var x{1..3} integer >= -3 <= 2;

maximize o: 0*x[1] + 0*x[2] + 4*x[3];

s.t. c1: max(x[2], x[1]) - min(x[2], x[1]) < -2;
s.t. c2: max(x[2], x[1]) - min(x[1], x[2]) >= -2;
s.t. c3: (x[2] <= 2) or (x[2] < 0) or (x[1] > -2);
