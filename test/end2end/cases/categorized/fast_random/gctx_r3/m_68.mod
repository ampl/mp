var x{1..5} integer >= 0 <= 3;

minimize o: 3*x[1] + 6*x[2] + 5*x[3] + 6*x[4] + 6*x[5];

s.t. c1:
   (if max(x[4], x[5]) != 2 then 1 else 0)
   + (if max(x[4], x[5]) != 2 then 1 else 0)
   + (if min(x[5], x[1]) > 1 then 1 else 0)
   >= 2;
s.t. c2: (min(x[5], x[1]) > 1) <==> (max(x[4], x[5]) != 2);
s.t. c3:
   (max(x[4], x[5]) != 2) ==> ((max(x[4], x[5]) != 2) or (min(x[5], x[1]) > 1));
s.t. c4:
   (if min(x[5], x[1]) > 1 then x[1] else 0)
   + (if min(x[5], x[1]) > 1 then x[1] else 0)
   > 0;
