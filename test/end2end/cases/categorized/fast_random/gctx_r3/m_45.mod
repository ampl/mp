var x{1..3} integer >= -1 <= 1;

maximize o: 4*x[1] - 5*x[2] - 3*x[3];

s.t. c1: (if min(x[1], x[3]) != -1 then x[2] else x[1]) == 1;
s.t. c2: not (min(x[1], x[3]) >= -1);
s.t. c3: not (max(x[1], x[2]) == -1);
s.t. c4:
   (if max(x[1], x[2]) == -1 then 1 else 0)
   + (if max(x[1], x[2]) == -1 then 1 else 0)
   + (if min(x[1], x[3]) >= -1 then 1 else 0)
   >= 3;
