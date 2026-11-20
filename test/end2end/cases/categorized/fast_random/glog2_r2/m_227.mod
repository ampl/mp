var x{1..6} integer >= 0 <= 2;

maximize o: -999*x[1] - 2000*x[2] + 2000*x[3] - 999*x[4] - 999*x[5] + 17*x[6];

s.t. c1: alldiff(x[2], x[1], x[4], x[6]);
s.t. c2: exists {i in 1..6} x[i] == 1;
s.t. c3: ((x[5] > 2) ==> (x[6] == 1)) and ((x[6] == 1) ==> (x[6] >= 2));
s.t. c4: (if x[3] < 1 then x[2] else x[3]) == 0;
