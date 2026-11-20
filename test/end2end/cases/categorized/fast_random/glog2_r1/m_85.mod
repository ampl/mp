var x{1..5} integer >= -2 <= 2;
minimize o: (0)*x[1] + (3)*x[2] + (-4)*x[3] + (4)*x[4] + (4)*x[5];
s.t. c1: ((if x[1] == 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
s.t. c2: alldiff(x[3], x[5], x[1]);
s.t. c3: (if x[5] != -2 then x[1] else x[4]) < 1;
s.t. c4: (x[4] > -1) and (x[5] < 1) ==> (x[1] >= -1);
s.t. c5: (x[5] != 0) ==> (x[3] != 0) else (x[3] == 2);
