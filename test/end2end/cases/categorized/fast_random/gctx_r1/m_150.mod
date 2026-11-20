var x{1..5} integer >= 0 <= 3;

maximize o: 2*x[1] + 4*x[2] + 3*x[3] + 0*x[4] + x[5];

s.t. c1:
   (if x[5] != x[4] then 1 else 0) + (if x[1] + x[2] != 0 then 1 else 0)
   + (if x[1] < 1 then 1 else 0)
   >= 2;
s.t. c2: (x[1] + x[2] != 0) ==> ((x[1] < 1) or (x[1] < 1));
