var x{1..4} integer >= 0 <= 3;

minimize o: 2*x[1] + 6*x[2] - 5*x[3] - 4*x[4];

s.t. c1: (max(x[1], x[2]) != 3) ==> ((max(x[1], x[2]) != 3) or (x[3] == x[3]));
s.t. c2: not (x[3] == x[3]);
s.t. c3: (if x[3] == x[3] then x[2] else x[3]) < 1;
s.t. c4:
   (if max(x[1], x[4]) <= 1 then 1 else 0) + (if x[3] == x[3] then 1 else 0)
   + (if x[3] == x[3] then 1 else 0)
   >= 1;
