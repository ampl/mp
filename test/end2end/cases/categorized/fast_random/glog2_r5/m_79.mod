var x{1..3} integer >= -2 <= 0;
maximize o: (3)*x[1] + (-5)*x[2] + (4)*x[3];
s.t. c1: (if x[3] < 0 then x[1] else x[2]) != 0;
s.t. c2: max(x[1], x[2]) - min(x[2], x[3]) > -1;
s.t. c3: (x[2] != -2) or (x[2] > 0) or (x[1] < -1);
s.t. c4: (if x[1] >= 0 then x[2] else x[3]) < -2;
s.t. c5: (if x[3] < 0 then x[1] else x[2]) != 0;
s.t. c6: (if x[1] >= 0 then x[2] else x[3]) < -2;
