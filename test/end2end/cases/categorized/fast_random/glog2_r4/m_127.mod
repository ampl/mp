var x{1..6} integer >= -2 <= 2;

maximize o: -2000*x[1] - 999*x[2] + 34*x[3] - 333*x[4] + 3000*x[5] - 666*x[6];

s.t. c1: (x[5] >= -2) ==> (x[5] == -1) else (x[6] >= 1);
s.t. c2: max(x[3], x[6]) - min(x[4], x[4]) == 1;
s.t. c3: (if x[4] <= -2 then x[5] else x[2]) == 1;
