var x{1..5} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (4)*x[2] + (5)*x[3] + (-4)*x[4] + (1)*x[5] + (4)*(if abs(x[5] - x[2]) == 1 then 1 else 0);
s.t. c1: (abs(x[5] - x[2]) == 1) and ((x[2] + x[5] > 0) or not (abs(x[5] - x[2]) == 1));
s.t. c2: not (abs(x[5] - x[2]) == 1);
s.t. c3: (abs(x[5] - x[2]) == 1) ==> (x[2] + x[5] > 0);
s.t. c4: (abs(x[5] - x[2]) == 1) ==> (abs(x[5] - x[2]) == 1) else (x[2] + x[5] > 0);
