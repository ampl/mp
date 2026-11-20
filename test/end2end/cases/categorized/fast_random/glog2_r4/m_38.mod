var x{1..5} integer >= 0 <= 1;
minimize o: (-4)*x[1] + (5)*x[2] + (-3)*x[3] + (5)*x[4] + (-3)*x[5];
s.t. c1: abs(x[4] - x[3]) != 1;
s.t. c2: (x[1] < 0) <==> (x[1] + x[2] <= 1);
s.t. c3: ((if x[3] >= 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
