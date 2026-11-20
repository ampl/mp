var x{1..4} integer >= -1 <= 1;

minimize o: 5*x[1] - 3*x[2] + 0*x[3] + x[4];

s.t. c1:
   (if x[1] != x[3] then x[4] else 0) + (if x[1] != x[3] then x[2] else 0) < -1;
s.t. c2: (x[3] <= x[3]) ==> (x[4] == 0);
s.t. c3: not ((x[3] <= x[3]) and (x[4] == 0));
s.t. c4: not (x[1] != x[3]);
