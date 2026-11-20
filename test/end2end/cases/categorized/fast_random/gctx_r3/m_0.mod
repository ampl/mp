var x{1..5} integer >= 0 <= 2;

maximize o: -4*x[1] + 4*x[2] + 5*x[3] + 2*x[4] + 4*x[5];

s.t. c1: (x[4] + x[1] > 2) ==> (x[1] > 2);
s.t. c2: not (x[4] + x[1] > 2);
s.t. c3: (if x[4] + x[1] > 2 then x[1] else x[5]) > 0;
s.t. c4:
   (if x[4] + x[1] > 2 then x[4] else 0) + (if x[1] > 2 then x[3] else 0) > 2;
