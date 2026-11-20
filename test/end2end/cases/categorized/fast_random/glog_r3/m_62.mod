var x{1..3} integer >= 0 <= 4;

minimize o: -3*x[1] - 5*x[2] + 3*x[3];

s.t. c1: alldiff(x[1], x[3], x[2]);
s.t. c2: (x[1] != 1) ==> (x[1] < 2);
s.t. c3: alldiff(x[3], x[1]);
s.t. c4: (if x[2] <= 3 then x[2] else x[1]) > 4;
