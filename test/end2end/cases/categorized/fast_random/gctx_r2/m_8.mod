var x{1..4} integer >= -1 <= 2;

maximize o: -5*x[1] - x[2] + 2*x[3] - x[4];

s.t. c1: (x[1] + x[3] >= 2) ==> (abs(x[4] - x[4]) != 1) else (x[2] > -1);
s.t. c2: (abs(x[4] - x[4]) != 1) ==> ((x[2] > -1) or (min(x[1], x[2]) >= 2));
s.t. c3: (abs(x[4] - x[4]) != 1) ==> (x[2] > -1) else (x[1] + x[3] >= 2);
s.t. c4:
   (if min(x[1], x[2]) >= 2 then x[1] else 0)
   + (if x[1] + x[3] >= 2 then x[2] else 0)
   > 0;
