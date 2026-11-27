var x{1..5} integer >= 0 <= 3;

minimize o: 6*x[1] - 2*x[2] - 5*x[3] + 4*x[4] + x[5];

s.t. c1: (x[1] >= 3) <==> (x[1] >= 3);
s.t. c2:
   (if max(x[3], x[1]) != 2 then 1 else 0)
   + (if max(x[3], x[1]) != 2 then 1 else 0)
   + (if max(x[3], x[1]) != 2 then 1 else 0)
   >= 1;
s.t. c3: not (x[1] >= 3);
s.t. c4: (if x[1] >= 3 then x[3] else x[4]) == 1;
