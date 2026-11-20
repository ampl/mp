var x{1..5} integer >= -1 <= 2;
minimize o: (4)*x[1] + (4)*x[2] + (6)*x[3] + (-4)*x[4] + (-5)*x[5] + (1)*(if x[2] < x[5] then 1 else 0);
s.t. c1: not ((x[1] >= x[5]) and (x[1] >= x[5]));
s.t. c2: (x[2] < x[5]) and ((x[2] < x[5]) or not (x[2] < x[5]));
s.t. c3: (x[2] < x[5]) ==> (x[1] >= x[5]) else (x[2] < x[5]);
