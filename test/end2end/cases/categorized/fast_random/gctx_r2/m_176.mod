var x{1..3} integer >= 0 <= 2;

minimize o: 2*x[1] - 5*x[2] + 3*x[3];

s.t. c1: (x[2] != x[2]) ==> (x[2] != x[2]) else (x[2] != x[2]);
s.t. c2: (x[2] != x[2]) or (x[2] != x[2]);
s.t. c3: (max(x[3], x[2]) == 0) ==> (max(x[3], x[2]) == 0) else (x[3] <= x[1]);
s.t. c4:
   (if max(x[3], x[2]) == 0 then 1 else 0) + (if x[2] != x[2] then 1 else 0)
   + (if max(x[3], x[2]) == 0 then 1 else 0)
   >= 2;
