var x{1..6} integer >= 0 <= 3;

minimize o: -333*x[1] + 17*x[2] - 333*x[3] + 17*x[4] + 17*x[5] - 2000*x[6];

s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) >= 3;
