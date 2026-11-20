var x{1..4} integer >= 0 <= 1;

minimize o: 2*x[1] + 6*x[2] + 6*x[3] + 3*x[4];

s.t. c1: max(x[1], x[2]) - min(x[3], x[3]) != 0;
s.t. c2: alldiff(x[2], x[3], x[1]);
s.t. c3: (x[3] < 1) or (x[3] != 1) or (x[1] > 1);
