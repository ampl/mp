var x{1..4} integer >= 0 <= 3;

maximize o: 0*x[1] + 0*x[2] + x[3] + x[4];

s.t. c1: (abs(x[1] - x[2]) >= 3) ==> (max(x[4], x[3]) >= 0);
s.t. c2: (if x[3] < 1 then x[2] else x[4]) == 2;
s.t. c3:
   (abs(x[1] - x[2]) >= 3)
   ==> ((abs(x[1] - x[2]) >= 3) or (max(x[4], x[3]) >= 0));
s.t. c4: not (max(x[4], x[3]) >= 0);
