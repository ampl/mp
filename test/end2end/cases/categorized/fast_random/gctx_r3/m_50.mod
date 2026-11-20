var x{1..3} integer >= 0 <= 2;

maximize o: -3*x[1] - 3*x[2] + 4*x[3];

s.t. c1:
   (if x[1] != 0 then 1 else 0) + (if x[1] > 1 then 1 else 0)
   + (if x[1] != 0 then 1 else 0)
   >= 1;
s.t. c2: (abs(x[1] - x[2]) >= 1) <==> (abs(x[1] - x[2]) >= 1);
s.t. c3: not (x[1] > 1);
