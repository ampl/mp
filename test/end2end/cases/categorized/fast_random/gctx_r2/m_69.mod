var x{1..3} integer >= 0 <= 3;

maximize o: -3*x[1] + x[2] + 2*x[3];

s.t. c1:
   (if max(x[2], x[2]) > 0 then 1 else 0)
   + (if max(x[2], x[2]) > 0 then 1 else 0)
   + (if x[1] + x[3] != 0 then 1 else 0)
   >= 1;
s.t. c2: not (max(x[2], x[2]) > 0);
s.t. c3: (if x[1] + x[3] != 0 then x[3] else x[2]) > 3;
s.t. c4: (if max(x[2], x[2]) > 0 then x[2] else x[3]) == 2;
