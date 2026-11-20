var x{1..3} integer >= -1 <= 1;

maximize o: x[1] - x[2] - x[3];

s.t. c1:
   (abs(x[1] - x[1]) <= 0) ==> ((abs(x[1] - x[1]) <= 0) or (x[2] <= x[2]));
s.t. c2: (x[2] <= x[2]) ==> (abs(x[1] - x[1]) <= 0) else (x[2] <= x[2]);
s.t. c3: (abs(x[1] - x[1]) <= 0) ==> (x[2] <= x[2]);
