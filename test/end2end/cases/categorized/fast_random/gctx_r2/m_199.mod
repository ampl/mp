var x{1..5} integer >= 0 <= 3;
maximize o: (0)*x[1] + (2)*x[2] + (-5)*x[3] + (2)*x[4] + (3)*x[5];
s.t. c1: (if x[3] != 3 then 1 else 0) + (if x[3] != 3 then 1 else 0) + (if x[2] > 3 then 1 else 0) >= 3;
s.t. c2: (x[3] != 3) ==> (x[3] != 3);
s.t. c3: (if x[2] < x[5] then 1 else 0) + (if x[3] != 3 then 1 else 0) + (if x[2] < x[5] then 1 else 0) >= 1;
s.t. c4: (x[2] > 3) <==> (x[3] != 3);
