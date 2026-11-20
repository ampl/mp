var x{1..5} integer >= 0 <= 3;

minimize o: -2*x[1] + 2*x[2] - 5*x[3] + 2*x[4] - 5*x[5];

s.t. c1:
   (if x[5] <= 2 then x[2] else 0) + (if x[4] + x[3] >= 1 then x[3] else 0) > 2;
s.t. c2: (x[5] <= 2) ==> (x[4] + x[3] >= 1);
