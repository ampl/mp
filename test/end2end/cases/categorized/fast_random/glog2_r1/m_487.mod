var x{1..3} integer >= 0 <= 2;

minimize o: -3000*x[1] - 333*x[2] - 333*x[3];

s.t. c1: max(x[3], x[2]) - min(x[2], x[1]) == 0;
