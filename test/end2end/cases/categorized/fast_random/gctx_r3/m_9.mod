var x{1..5} integer >= -1 <= 1;

minimize o:
   -3*x[1] - 3*x[2] - x[3] + 6*x[4] + x[5]
   + 6*(if max(x[3], x[1]) == 1 then 1 else 0);

s.t. c1:
   (if x[5] >= -1 then x[2] else 0) + (if x[5] >= -1 then x[5] else 0) >= 1;
s.t. c2: (x[5] >= -1) ==> (max(x[3], x[1]) == 1) else (x[5] >= -1);
s.t. c3: (x[5] >= -1) ==> (max(x[3], x[1]) == 1) else (x[5] >= -1);
s.t. c4: (max(x[3], x[1]) == 1) or (max(x[3], x[1]) == 1);
