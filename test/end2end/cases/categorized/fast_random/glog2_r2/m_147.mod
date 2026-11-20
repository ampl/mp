var x{1..5} integer >= 0 <= 3;

minimize o: 17*x[1] - 3000*x[2] - 999*x[3] + 34*x[4] - 1000*x[5];

s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) != 1;
