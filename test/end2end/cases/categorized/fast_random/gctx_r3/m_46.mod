var x{1..3} integer >= 0 <= 2;

maximize o: 2*x[1] - 4*x[2] + 3*x[3];

s.t. c1:
   (abs(x[3] - x[2]) == 2) and ((x[3] == x[3]) or not (max(x[3], x[1]) > 1));
s.t. c2: (x[3] == x[3]) <==> (x[3] == x[3]);
s.t. c3:
   (min(x[3], x[1]) >= 2)
   ==> (min(x[3], x[1]) >= 2)
   else (abs(x[3] - x[2]) == 2);
