var x{1..6} integer >= 0 <= 2;

minimize o:
   -666*x[1] - 2000*x[2] - 333*x[3] + 2000*x[4] + 1000*x[5] + 2000*x[6];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 1;
