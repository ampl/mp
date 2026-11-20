var x{1..4} integer >= 0 <= 2;

maximize o: 3*x[1] + 4*x[2] + 3*x[3] - 5*x[4];

s.t. c1:
   (if x[4] >= x[1] then 1 else 0) + (if x[4] >= x[1] then 1 else 0)
   + (if x[4] >= x[1] then 1 else 0)
   >= 3;
s.t. c2: (if x[4] >= x[1] then x[4] else x[2]) == 0;
s.t. c3:
   (min(x[1], x[3]) > 1) ==> (min(x[2], x[4]) == 0) else (min(x[1], x[3]) > 1);
s.t. c4: (x[4] >= x[1]) <==> (min(x[2], x[4]) == 0);
