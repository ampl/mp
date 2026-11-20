var x{1..6} integer >= -2 <= 2;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5] + 0*x[6];

s.t. c1: alldiff(x[6], x[3]);
s.t. c2: (if x[4] >= 2 then x[5] else x[3]) != -1;
