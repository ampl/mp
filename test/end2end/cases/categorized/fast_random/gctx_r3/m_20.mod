var x{1..3} integer >= 0 <= 3;

maximize o: -2*x[1] + 3*x[2] + 3*x[3];

s.t. c1:
   (max(x[3], x[2]) < 2) and ((x[1] == x[2]) or not (max(x[3], x[2]) < 2));
s.t. c2: (if x[1] == x[2] then x[1] else x[2]) > 0;
s.t. c3: (max(x[3], x[2]) < 2) or (x[1] == x[2]);
s.t. c4: not (x[1] == x[2]);
