var x{1..6} integer >= 0 <= 4;

maximize o: -x[1] - 2*x[2] - 5*x[3] + 3*x[4] - 5*x[5] - 4*x[6];

s.t. c1: alldiff(x[3], x[5]);
s.t. c2: (if x[2] == 2 then x[2] else x[1]) >= 4;
s.t. c3: max(x[1], x[1]) - min(x[4], x[1]) > 2;
s.t. c4: (x[2] > 2) ==> (x[4] > 4);
