var x{1..5} integer >= 0 <= 4;

maximize o: 6*x[1] + 3*x[2] + 6*x[3] + x[4] - x[5];

s.t. c1: max(x[5], x[3]) - min(x[2], x[2]) == 1;
s.t. c2: (if x[3] != 1 then x[4] else x[2]) <= 1;
s.t. c3: max(x[1], x[5]) - min(x[1], x[2]) > 4;
s.t. c4: (x[4] != 2) or (x[4] >= 3) or (x[5] <= 0);
s.t. c5: max(x[5], x[3]) - min(x[2], x[2]) == 1;
s.t. c6: (x[4] != 2) or (x[4] >= 3) or (x[5] <= 0);
