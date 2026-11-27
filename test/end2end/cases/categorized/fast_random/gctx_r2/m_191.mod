var x{1..5} integer >= 0 <= 2;

maximize o: 6*x[1] - 4*x[2] - 5*x[3] - 3*x[4] + 2*x[5];

s.t. c1: (x[3] >= x[5]) <==> (max(x[4], x[1]) != 0);
s.t. c2:
   (if x[3] >= x[5] then x[2] else 0) + (if x[3] >= x[5] then x[1] else 0) > 2;
s.t. c3: not ((max(x[4], x[1]) != 0) and (x[3] >= x[5]));
s.t. c4: (x[1] + x[2] != 1) ==> (x[3] >= x[5]);
