var x{1..6} integer >= 0 <= 3;

minimize o: 3*x[1] - 4*x[2] - 5*x[3] - 3*x[4] + 0*x[5] + 2*x[6];

s.t. c1: not ((x[5] >= 3) and (x[4] > 2));
s.t. c2: (x[5] == 3) and (x[1] >= 2) ==> (x[3] >= 3);
s.t. c3: (if x[4] > 0 then x[1] else x[1]) >= 2;
s.t. c4: (x[5] <= 2) or (x[3] <= 2) or (x[4] != 3);
s.t. c5: atmost 3 {i in 1..6} (x[i] >= 0);
s.t. c6: alldiff(x[1], x[3]);
