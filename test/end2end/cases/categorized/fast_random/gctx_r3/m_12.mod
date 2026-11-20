var x{1..4} integer >= 0 <= 3;
maximize o: (6)*x[1] + (6)*x[2] + (-3)*x[3] + (2)*x[4];
s.t. c1: (x[3] + x[3] < 1) ==> ((x[3] + x[3] < 1) or (x[3] + x[3] < 1));
s.t. c2: (x[3] + x[3] < 1) ==> ((x[3] + x[3] < 1) or (max(x[4], x[4]) >= 1));
s.t. c3: (if x[3] + x[3] < 1 then 1 else 0) + (if x[3] + x[3] < 1 then 1 else 0) + (if x[3] + x[3] < 1 then 1 else 0) >= 1;
s.t. c4: (x[3] + x[3] < 1) or (max(x[4], x[4]) >= 1);
