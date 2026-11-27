var x{1..4} integer >= 0 <= 3;

maximize o:
   3*x[1] + 2*x[2] - 5*x[3] - 4*x[4]
   + 2*(if abs(x[2] - x[1]) < 1 then 1 else 0);

s.t. c1: (if x[2] >= 1 then x[2] else x[1]) < 0;
s.t. c2: not (x[2] >= 1);
s.t. c3: (x[2] >= 1) ==> ((x[4] == x[2]) or (x[4] == x[2]));
s.t. c4: not ((abs(x[2] - x[1]) < 1) and (abs(x[2] - x[1]) < 1));
