var x{1..5} integer >= 0 <= 3;
maximize o: (-4)*x[1] + (-2)*x[2] + (4)*x[3] + (-4)*x[4] + (1)*x[5];
s.t. c1: (abs(x[2] - x[2]) < 1) or (abs(x[2] - x[2]) < 1);
s.t. c2: (if abs(x[2] - x[2]) < 1 then 1 else 0) + (if x[4] + x[4] != 2 then 1 else 0) + (if abs(x[2] - x[2]) < 1 then 1 else 0) >= 1;
s.t. c3: (abs(x[2] - x[3]) != 0) ==> (x[4] + x[4] != 2);
