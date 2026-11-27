var x{1..5} integer >= -1 <= 1;

minimize o: 2*x[1] + 4*x[2] + 6*x[3] + 4*x[4] - 2*x[5];

s.t. c1:
   (if x[4] < 1 then 1 else 0) + (if x[4] < 1 then 1 else 0)
   + (if x[2] + x[1] <= -1 then 1 else 0)
   >= 1;
s.t. c2: (x[2] > -1) ==> (x[2] > -1) else (x[2] > -1);
s.t. c3:
   (if max(x[3], x[5]) != 1 then x[2] else 0) + (if x[2] > -1 then x[2] else 0)
   <= 1;
s.t. c4: (x[2] + x[1] <= -1) or (x[2] > -1);
