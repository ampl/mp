var x{1..5} integer >= -1 <= 2;
minimize o: (2)*x[1] + (-4)*x[2] + (5)*x[3] + (1)*x[4] + (-3)*x[5];
s.t. c1: (if abs(x[3] - x[1]) > 2 then 1 else 0) + (if abs(x[3] - x[1]) > 2 then 1 else 0) + (if abs(x[3] - x[1]) > 2 then 1 else 0) >= 1;
s.t. c2: (x[1] < x[5]) ==> (abs(x[3] - x[1]) > 2) else (x[5] > 0);
