var x{1..5} integer >= 0 <= 3;

maximize o: -x[1] + 2*x[2] + x[3] - x[4] - 5*x[5];

s.t. c1: (x[4] > 3) and ((x[4] > 3) or not (x[4] > 3));
s.t. c2: (x[4] > 3) ==> (x[4] > 3);
s.t. c3:
   (if x[2] == x[4] then 1 else 0) + (if x[2] == x[4] then 1 else 0)
   + (if x[4] > 3 then 1 else 0)
   >= 2;
s.t. c4: (x[2] == x[4]) or (x[2] == x[4]);
