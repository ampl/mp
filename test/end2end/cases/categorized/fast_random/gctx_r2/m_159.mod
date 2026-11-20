var x{1..4} integer >= -1 <= 1;
maximize o: (0)*x[1] + (1)*x[2] + (3)*x[3] + (-1)*x[4];
s.t. c1: (if abs(x[2] - x[2]) < 1 then x[3] else 0) + (if abs(x[2] - x[2]) < 1 then x[4] else 0) > 1;
s.t. c2: (if x[1] + x[2] == 0 then 1 else 0) + (if x[2] != 1 then 1 else 0) + (if abs(x[2] - x[2]) < 1 then 1 else 0) >= 3;
s.t. c3: (x[2] != 1) ==> (abs(x[2] - x[2]) < 1) else (x[2] != 1);
