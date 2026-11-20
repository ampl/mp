var x{1..4} integer >= 0 <= 3;

maximize o: 5*x[1] - 5*x[2] - 3*x[3] + x[4] + 3*(if x[3] < 2 then 1 else 0);

s.t. c1: (x[3] < 2) <==> (x[3] < 2);
s.t. c2: (x[3] < 2) ==> ((x[3] < 2) or (min(x[2], x[1]) <= 0));
s.t. c3:
   (if min(x[2], x[1]) <= 0 then x[3] else 0) + (if x[3] < 2 then x[2] else 0)
   >= 0;
