var x{1..5} integer >= 0 <= 3;
minimize o: (4)*x[1] + (5)*x[2] + (-2)*x[3] + (5)*x[4] + (-5)*x[5];
s.t. c1: ((if x[2] > 1 then 1 else 0) + (if x[4] >= 1 then 1 else 0)) >= 1;
s.t. c2: ((if x[1] >= 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
s.t. c3: (x[1] < 2) <==> (x[3] + x[1] >= 2);
s.t. c4: (x[3] <= 1) or (x[3] >= 1) or (x[2] <= 2);
s.t. c5: x[1] == 2;
