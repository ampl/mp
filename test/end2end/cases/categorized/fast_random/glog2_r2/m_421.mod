var x{1..3} integer >= -2 <= 1;

minimize o: x[1] + 5*x[2] + 5*x[3];

s.t. c1: numberof -3 in (x[1], x[2], x[3]) < 3;
