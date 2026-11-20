var x{1..3} integer >= 0 <= 3;
maximize o: (-5)*x[1] + (-1)*x[2] + (-5)*x[3];
s.t. c1: (x[2] + x[3] < 1) ==> (abs(x[1] - x[2]) >= 3) else (x[1] <= x[2]);
s.t. c2: (if x[2] + x[3] < 1 then 1 else 0) + (if abs(x[1] - x[2]) >= 3 then 1 else 0) + (if x[2] + x[3] < 1 then 1 else 0) >= 3;
