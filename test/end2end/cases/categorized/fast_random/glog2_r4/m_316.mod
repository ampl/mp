var x{1..3} integer >= 0 <= 2;

maximize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: (if x[1] > 1 then x[1] else x[2]) == 1;
s.t. c2: alldiff(x[3], x[1], x[2]);
