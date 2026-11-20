var x{1..5} integer >= -1 <= 2;

maximize o:
   -4*x[1] + 2*x[2] - 5*x[3] - 4*x[4] + 3*x[5] + 4*(if x[5] < 1 then 1 else 0);

s.t. c1:
   (if x[4] != 1 then 1 else 0) + (if x[4] != 1 then 1 else 0)
   + (if x[5] < 1 then 1 else 0)
   >= 2;
s.t. c2: (x[5] < 2) ==> (x[4] != 1);
s.t. c3: (x[5] < 1) or (x[5] < 2);
s.t. c4: (x[5] < 1) ==> ((x[5] < 1) or (x[5] < 2));
