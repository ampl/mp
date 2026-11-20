var x{1..6} integer >= -2 <= 1;

maximize o: 2000*x[1] - 666*x[2] + 51*x[3] + 34*x[4] - 333*x[5] - 1000*x[6];

s.t. c1: max(x[6], x[6]) - min(x[2], x[5]) < 1;
