var x{1..3} integer >= -1 <= 1;

maximize o: -3*x[1] - 2*x[2] + 4*x[3];

s.t. c1:
   (x[1] + x[3] >= 1) and ((x[1] + x[3] >= 1) or not (abs(x[2] - x[2]) >= 2));
s.t. c2: (if abs(x[2] - x[2]) >= 2 then x[3] else x[1]) == 0;
s.t. c3: (x[1] + x[3] >= 1) ==> (x[1] + x[3] >= 1) else (abs(x[2] - x[2]) >= 2);
