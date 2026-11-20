var x{1..4} integer >= 0 <= 3;
minimize o: (-4)*x[1] + (-5)*x[2] + (5)*x[3] + (-5)*x[4];
s.t. c1: (if x[2] >= 3 then x[4] else x[3]) < 3;
s.t. c2: not ((x[1] != 3) and (x[4] <= 0));
s.t. c3: ((x[1] <= 2) ==> (x[4] < 3)) and ((x[4] < 3) ==> (x[1] >= 2));
s.t. c4: not ((x[2] <= 3) and (x[2] > 1));
