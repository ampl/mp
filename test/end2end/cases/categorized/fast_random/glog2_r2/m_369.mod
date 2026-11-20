var x{1..3} integer >= 0 <= 2;
minimize o: (3)*x[1] + (5)*x[2] + (4)*x[3];
s.t. c1: ((if x[3] <= 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c2: (if x[1] > 0 then x[3] else x[1]) > 2;
s.t. c3: (x[3] >= 0) <==> (x[1] + x[3] == 2);
s.t. c4: abs(x[1] - x[3]) <= 2;
s.t. c5: ((if x[3] <= 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c6: abs(x[1] - x[3]) <= 2;
