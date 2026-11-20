var x{1..3} integer >= 0 <= 2;

maximize o: x[1] + 4*x[2] - 4*x[3];

s.t. c1: alldiff(x[2], x[3], x[1]);
s.t. c2: (if x[3] >= 2 then x[3] else x[2]) != 2;
s.t. c3: exists {i in 1..3} x[i] == 0;
s.t. c4: (x[3] > 0) and (x[1] < 2) ==> (x[3] >= 0);
s.t. c5: max(x[1], x[2]) - min(x[1], x[2]) > 0;
