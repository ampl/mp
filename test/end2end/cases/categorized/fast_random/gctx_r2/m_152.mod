var x{1..4} integer >= 0 <= 2;

minimize o: -x[1] + 2*x[2] + 2*x[3] + 4*x[4];

s.t. c1:
   (max(x[4], x[1]) < 1) ==> ((max(x[3], x[1]) >= 0) or (max(x[4], x[1]) < 1));
s.t. c2:
   (max(x[3], x[1]) >= 0) ==> (x[3] + x[1] == 1) else (max(x[4], x[1]) < 1);
s.t. c3:
   (if max(x[3], x[1]) >= 0 then 1 else 0)
   + (if max(x[3], x[1]) >= 0 then 1 else 0)
   + (if x[3] + x[1] == 1 then 1 else 0)
   >= 1;
s.t. c4: not (max(x[3], x[1]) >= 0);
