var x{1..5} integer >= 0 <= 2;

maximize o: 0*x[1] - 5*x[2] - x[3] - 2*x[4] + 4*x[5];

s.t. c1: (x[3] + x[5] >= 2) ==> (x[4] < x[1]);
s.t. c2:
   (if x[5] < x[1] then x[5] else 0) + (if x[1] + x[1] >= 2 then x[4] else 0)
   >= 1;
s.t. c3: (x[4] < x[1]) ==> ((x[1] + x[1] >= 2) or (x[4] < x[1]));
s.t. c4: (x[5] < x[1]) <==> (x[3] + x[5] >= 2);
