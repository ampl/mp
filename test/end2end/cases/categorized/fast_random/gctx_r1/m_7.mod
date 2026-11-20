var x{1..3} integer >= -1 <= 2;
minimize o: (4)*x[1] + (5)*x[2] + (-4)*x[3];
s.t. c1: (if x[3] < 1 then x[2] else 0) + (if x[3] < 1 then x[3] else 0) > 2;
s.t. c2: (x[3] < 1) ==> (x[3] < 1);
s.t. c3: (x[1] < x[3]) ==> (x[1] < x[3]) else (x[1] < x[3]);
s.t. c4: (if x[3] < 1 then 1 else 0) + (if x[3] < 1 then 1 else 0) + (if x[3] < 1 then 1 else 0) >= 1;
