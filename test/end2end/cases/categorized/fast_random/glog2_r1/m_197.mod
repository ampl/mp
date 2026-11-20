var x{1..3} integer >= 0 <= 2;

minimize o: -666*x[1] - 2000*x[2] + 3000*x[3];

s.t. c1: not ((x[1] > 2) and (x[2] <= 2));
