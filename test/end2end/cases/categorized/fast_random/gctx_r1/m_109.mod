var x{1..5} integer >= -1 <= 2;

maximize o: 5*x[1] + 3*x[2] + 6*x[3] + 4*x[4] + x[5];

s.t. c1:
   (min(x[4], x[3]) > -1)
   ==> ((min(x[2], x[1]) <= 0) or (min(x[2], x[1]) <= 0));
s.t. c2: (min(x[4], x[3]) > -1) ==> (min(x[2], x[1]) <= 0);
