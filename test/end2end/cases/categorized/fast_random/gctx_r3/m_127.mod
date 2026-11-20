var x{1..5} integer >= -1 <= 2;

minimize o:
   -5*x[1] - 5*x[2] + 0*x[3] - 5*x[4] - 5*x[5]
   + 4*(if min(x[2], x[1]) != -1 then 1 else 0);

s.t. c1:
   (if x[5] + x[4] < -1 then x[3] else 0)
   + (if x[5] + x[4] < -1 then x[3] else 0)
   >= 0;
s.t. c2: (x[5] + x[4] < -1) <==> (min(x[2], x[1]) != -1);
s.t. c3: (x[5] + x[4] < -1) ==> ((x[5] + x[4] < -1) or (min(x[2], x[1]) != -1));
s.t. c4:
   (min(x[2], x[1]) != -1) ==> ((x[5] + x[4] < -1) or (min(x[2], x[1]) != -1));
