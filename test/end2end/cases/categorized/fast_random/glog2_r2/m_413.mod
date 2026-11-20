var x{1..3} integer >= -2 <= 0;

maximize o: -3*x[1] - x[2] + 3*x[3];

s.t. c1: max(x[3], x[3]) - min(x[1], x[2]) >= -2;
s.t. c2: x[3] == -1;
