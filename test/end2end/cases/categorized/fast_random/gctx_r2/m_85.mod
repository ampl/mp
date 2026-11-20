var x{1..5} integer >= -1 <= 2;

maximize o: -2*x[1] + 0*x[2] - 3*x[3] + 5*x[4] + 5*x[5];

s.t. c1:
   (if x[2] + x[1] == 0 then 1 else 0) + (if x[4] != 2 then 1 else 0)
   + (if x[2] + x[1] == 0 then 1 else 0)
   >= 1;
s.t. c2: (x[2] + x[1] == 0) and ((x[4] != 2) or not (x[4] != 2));
