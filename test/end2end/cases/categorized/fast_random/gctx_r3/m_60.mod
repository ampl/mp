var x{1..3} integer >= -1 <= 2;
minimize o: (-1)*x[1] + (5)*x[2] + (1)*x[3];
s.t. c1: (x[3] < x[3]) and ((x[3] == x[1]) or not (x[2] > 0));
s.t. c2: (x[3] == x[1]) ==> (x[3] < x[3]);
s.t. c3: (if x[3] < x[3] then x[3] else x[2]) >= 1;
s.t. c4: (x[3] < x[3]) or (x[3] == x[1]);
