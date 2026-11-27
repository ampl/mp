var x{1..6} integer >= 0 <= 2;

maximize o: -4*x[1] + 3*x[2] - x[3] - 4*x[4] - 3*x[5] - 3*x[6];

s.t. c1: (x[5] > 2) ==> (x[1] < 2) else (x[6] == 1);
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) == 3;
s.t. c3: (if x[1] < 0 then x[1] else x[3]) < 2;
s.t. c4: not ((x[2] > 0) and (x[1] != 0));
