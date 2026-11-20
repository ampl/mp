var x{1..5} integer >= -1 <= 2;

maximize o: 6*x[1] - 5*x[2] + 2*x[3] + x[4] + 3*x[5];

s.t. c1: (min(x[4], x[5]) < 1) <==> (min(x[4], x[5]) < 1);
s.t. c2:
   (if min(x[4], x[5]) < 1 then 1 else 0) + (if x[1] >= x[4] then 1 else 0)
   + (if x[3] != 1 then 1 else 0)
   >= 3;
s.t. c3: (x[1] >= x[4]) ==> (min(x[4], x[5]) < 1) else (x[1] >= x[4]);
s.t. c4: (x[1] >= x[4]) and ((x[3] != 1) or not (x[3] != 1));
