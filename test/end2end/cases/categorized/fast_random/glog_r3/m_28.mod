var x{1..5} integer >= 0 <= 3;

maximize o: 0*x[1] + 5*x[2] + 3*x[3] + 5*x[4] + 6*x[5];

s.t. c1: (if x[5] != 3 then x[2] else x[5]) >= 1;
s.t. c2: alldiff(x[1], x[3]);
s.t. c3: max(x[4], x[2]) - min(x[5], x[5]) >= 0;
s.t. c4: max(x[2], x[1]) - min(x[2], x[4]) >= 1;
