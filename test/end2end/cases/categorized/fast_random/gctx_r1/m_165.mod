var x{1..4} integer >= 0 <= 2;

minimize o: 5*x[1] + 5*x[2] + 5*x[3] - 4*x[4];

s.t. c1:
   (if x[2] < x[2] then x[4] else 0) + (if x[2] < x[2] then x[1] else 0) >= 2;
s.t. c2: (max(x[1], x[3]) >= 2) ==> (x[2] < x[2]) else (x[2] < x[2]);
s.t. c3: (max(x[1], x[3]) >= 2) ==> (max(x[1], x[3]) >= 2);
