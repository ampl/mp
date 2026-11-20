var x{1..5} integer >= 0 <= 2;

maximize o:
   4*x[1] - 2*x[2] + 6*x[3] - x[4] - 3*x[5] - (if x[4] != 0 then 1 else 0);

s.t. c1: (if x[1] == x[1] then x[3] else x[1]) <= 2;
s.t. c2: (if x[4] != 0 then x[2] else x[3]) <= 2;
s.t. c3: (x[1] == x[1]) ==> (abs(x[4] - x[2]) <= 2);
