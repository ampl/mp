var x{1..3} integer >= -2 <= 0;

maximize o: x[1] + 5*x[2] + 4*x[3];

s.t. c1: ((if x[1] != -2 then 1 else 0) + (if x[1] >= -2 then 1 else 0)) >= 1;
s.t. c2: ((x[1] != 0) ==> (x[2] <= -1)) and ((x[2] <= -1) ==> (x[2] >= 0));
s.t. c3: (x[3] < -1) <==> (x[3] + x[2] > -2);
