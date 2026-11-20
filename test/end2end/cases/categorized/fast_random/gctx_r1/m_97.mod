var x{1..5} integer >= -1 <= 1;

maximize o: 2*x[1] + 2*x[2] - 3*x[3] + 0*x[4] - x[5];

s.t. c1: not (abs(x[2] - x[2]) == 1);
s.t. c2:
   (if abs(x[2] - x[2]) == 1 then x[2] else 0) + (if x[1] < 0 then x[4] else 0)
   > -1;
s.t. c3: not (x[1] < 0);
s.t. c4: (x[1] < 0) ==> (abs(x[2] - x[2]) == 1) else (abs(x[2] - x[2]) == 1);
