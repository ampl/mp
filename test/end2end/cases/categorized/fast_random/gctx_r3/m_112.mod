var x{1..4} integer >= 0 <= 3;
maximize o: (6)*x[1] + (4)*x[2] + (0)*x[3] + (-5)*x[4] + (0)*(if abs(x[2] - x[1]) < 2 then 1 else 0);
s.t. c1: (abs(x[2] - x[1]) < 2) and ((abs(x[2] - x[1]) < 2) or not (x[4] <= 3));
s.t. c2: (x[1] + x[2] != 2) ==> ((x[4] <= 3) or (abs(x[2] - x[1]) < 2));
s.t. c3: (x[4] <= 3) ==> ((x[4] <= 3) or (abs(x[2] - x[1]) < 2));
