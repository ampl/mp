var x{1..3} integer >= 0 <= 3;

minimize o: -x[1] + 3*x[2] + 6*x[3];

s.t. c1: ((if x[1] != 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
s.t. c2: (if x[1] > 1 then x[3] else x[3]) != 0;
s.t. c3: (x[2] == 2) <==> (x[2] + x[3] >= 3);
s.t. c4: ((x[3] >= 2) ==> (x[3] < 1)) and ((x[3] < 1) ==> (x[1] >= 2));
s.t. c5: x[2] == 0;
s.t. c6: x[3] == 0;
