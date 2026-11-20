var x{1..3} integer >= 0 <= 2;

maximize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: max(x[3], x[3]) - min(x[1], x[3]) >= 0;
