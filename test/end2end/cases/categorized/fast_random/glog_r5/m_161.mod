var x{1..3} integer >= 0 <= 3;

maximize o: 3*x[1] + x[2] + 4*x[3];

s.t. c1: alldiff(x[1], x[3], x[2]);
s.t. c2: (if x[2] != 2 then x[1] else x[3]) < 2;
