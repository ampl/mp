var x{1..4} integer >= 0 <= 3;

maximize o: 2*x[1] + 6*x[2] + 0*x[3] - 2*x[4];

s.t. c1:
   (if x[4] < 0 then 1 else 0) + (if x[3] >= 3 then 1 else 0)
   + (if x[3] >= 3 then 1 else 0)
   >= 1;
s.t. c2: (x[1] != 3) ==> (min(x[2], x[3]) <= 2);
s.t. c3: (x[3] >= 3) <==> (x[4] < 0);
s.t. c4: (if x[3] >= 3 then x[1] else x[1]) == 3;
