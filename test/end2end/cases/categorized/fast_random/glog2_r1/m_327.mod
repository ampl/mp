var x{1..6} integer >= 0 <= 2;

maximize o:
   -666*x[1] - 333*x[2] + 2000*x[3] - 1000*x[4] + 3000*x[5] + 3000*x[6];

s.t. c1: alldiff(x[5], x[1], x[4], x[2]);
