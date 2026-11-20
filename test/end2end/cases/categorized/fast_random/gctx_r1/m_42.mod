var x{1..5} integer >= 0 <= 3;

maximize o: 4*x[1] - 2*x[2] + 2*x[3] + 2*x[4] + x[5];

s.t. c1:
   (if min(x[1], x[4]) <= 3 then 1 else 0) + (if x[4] != x[1] then 1 else 0)
   + (if x[4] != x[1] then 1 else 0)
   >= 3;
s.t. c2: (min(x[1], x[4]) <= 3) <==> (min(x[1], x[4]) <= 3);
s.t. c3: (min(x[4], x[3]) != 0) ==> (x[4] != x[1]) else (x[4] != x[1]);
