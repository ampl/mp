var x{1..3} integer >= 0 <= 2;

minimize o: 3*x[1] - 4*x[2] + 3*x[3] + 2*(if x[2] > x[3] then 1 else 0);

s.t. c1: (x[2] >= x[2]) <==> (max(x[1], x[1]) >= 0);
s.t. c2: (x[2] > x[3]) ==> ((x[2] > x[3]) or (x[2] > x[3]));
s.t. c3:
   (if max(x[1], x[1]) >= 0 then x[3] else 0)
   + (if max(x[1], x[1]) >= 0 then x[1] else 0)
   < 0;
s.t. c4:
   (if x[3] + x[1] != 1 then x[1] else 0)
   + (if max(x[1], x[1]) >= 0 then x[1] else 0)
   <= 2;
