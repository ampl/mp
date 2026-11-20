var x{1..5} integer >= -1 <= 2;
minimize o: (-2)*x[1] + (4)*x[2] + (2)*x[3] + (-5)*x[4] + (-2)*x[5];
s.t. c1: not (x[2] == 1);
s.t. c2: (x[2] == 1) ==> (x[5] != 2) else (x[2] == 1);
s.t. c3: (if x[2] == 1 then x[4] else x[4]) <= -1;
s.t. c4: (x[5] != 2) and ((x[2] == 1) or not (x[2] == 1));
