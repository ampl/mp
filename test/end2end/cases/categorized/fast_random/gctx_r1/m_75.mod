var x{1..5} integer >= -1 <= 1;
minimize o: (4)*x[1] + (4)*x[2] + (-1)*x[3] + (-3)*x[4] + (-4)*x[5];
s.t. c1: (max(x[4], x[5]) < -1) and ((max(x[4], x[5]) < -1) or not (x[3] != x[3]));
s.t. c2: (x[3] != x[3]) ==> (x[3] != x[3]) else (x[3] != x[3]);
