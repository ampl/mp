var x{1..3} integer >= 0 <= 2;
minimize o: (-1)*x[1] + (3)*x[2] + (4)*x[3];
s.t. c1: (x[3] <= 1) ==> (x[3] < 2);
s.t. c2: (x[2] == 2) <==> (x[2] + x[1] >= 0);
s.t. c3: ((if x[1] < 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
s.t. c4: (if x[1] == 0 then x[3] else x[2]) < 2;
s.t. c5: (x[3] <= 1) ==> (x[3] < 2);
s.t. c6: (if x[1] == 0 then x[3] else x[2]) < 2;
