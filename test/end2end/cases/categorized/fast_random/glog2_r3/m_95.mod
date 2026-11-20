var x{1..6} integer >= 0 <= 4;

minimize o: -2*x[1] + 6*x[2] + 0*x[3] - x[4] + 0*x[5] - 3*x[6];

s.t. c1: (if x[3] > 3 then x[4] else x[6]) != 3;
s.t. c2: not ((x[2] != 4) and (x[5] < 0));
s.t. c3: (x[4] == 1) or (x[6] == 3) or (x[1] == 0);
s.t. c4: max(x[6], x[3]) - min(x[4], x[1]) == 2;
s.t. c5: alldiff(x[1], x[4], x[6]);
s.t. c6: forall {i in 1..5} x[i] != x[i+1];
s.t. c7: (x[2] <= 1) or (x[2] != 2) or (x[6] <= 0);
