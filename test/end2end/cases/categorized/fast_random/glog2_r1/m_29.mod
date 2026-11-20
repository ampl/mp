var x{1..5} integer >= 0 <= 2;

minimize o: 2*x[1] + 0*x[2] - x[3] + 6*x[4] - 4*x[5];

s.t. c1: max(x[1], x[5]) - min(x[4], x[2]) >= 1;
s.t. c2: (if x[4] == 0 then x[3] else x[3]) != 0;
s.t. c3: (x[2] != 1) ==> (x[3] <= 0) else (x[5] >= 2);
s.t. c4: max(x[3], x[1]) - min(x[4], x[1]) == 0;
s.t. c5: max(x[1], x[5]) - min(x[4], x[2]) >= 1;
s.t. c6: max(x[3], x[1]) - min(x[4], x[1]) == 0;
