var x{1..3} integer >= -2 <= 2;
maximize o: (6)*x[1] + (-3)*x[2] + (-3)*x[3];
s.t. c1: (if x[3] < 2 then x[3] else x[1]) < -1;
s.t. c2: (x[3] == -2) ==> (x[3] > -2) else (x[1] <= -1);
s.t. c3: x[1] == 0;
