var x{1..5} integer >= 0 <= 3;

maximize o: 2*x[1] + 6*x[2] - x[3] + 3*x[4] - 2*x[5];

s.t. c1:
   (if x[4] > 0 then x[4] else 0) + (if abs(x[1] - x[1]) < 1 then x[5] else 0)
   < 2;
s.t. c2:
   (if x[4] > 0 then 1 else 0) + (if x[4] > 0 then 1 else 0)
   + (if abs(x[1] - x[1]) < 1 then 1 else 0)
   >= 1;
s.t. c3: (abs(x[1] - x[1]) < 1) ==> ((x[4] > 0) or (abs(x[1] - x[1]) < 1));
