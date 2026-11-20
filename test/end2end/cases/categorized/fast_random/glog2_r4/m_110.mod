var x{1..5} integer >= 0 <= 4;

maximize o: -x[1] + 2*x[2] + 4*x[3] + 6*x[4] + 6*x[5];

s.t. c1: max(x[2], x[2]) - min(x[1], x[3]) >= 2;
s.t. c2: alldiff(x[3], x[2]);
