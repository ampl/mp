var x{1..3} integer >= -2 <= 1;
minimize o: (2)*x[1] + (1)*x[2] + (6)*x[3];
s.t. c1: (x[3] >= -1) ==> (x[3] != 1) else (x[2] <= 0);
s.t. c2: exists {i in 1..3} x[i] == 0;
s.t. c3: max(x[2], x[1]) - min(x[2], x[2]) < -2;
s.t. c4: ((if x[3] > 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c5: max(x[3], x[2]) - min(x[2], x[1]) > -2;
s.t. c6: (x[1] < -2) and (x[3] >= -2) ==> (x[2] >= -2);
s.t. c7: max(x[3], x[3]) - min(x[3], x[1]) != -2;
