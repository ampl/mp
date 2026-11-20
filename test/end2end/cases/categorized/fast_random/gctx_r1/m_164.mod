var x{1..5} integer >= 0 <= 3;
minimize o: (-2)*x[1] + (0)*x[2] + (3)*x[3] + (0)*x[4] + (-5)*x[5] + (-4)*(if x[4] + x[3] == 2 then 1 else 0);
s.t. c1: not ((x[2] > x[4]) and (x[4] + x[3] == 2));
s.t. c2: (x[2] > x[4]) ==> ((x[2] > x[4]) or (x[2] > x[4]));
s.t. c3: (if x[4] + x[3] == 2 then x[2] else 0) + (if x[2] > x[4] then x[1] else 0) < 3;
