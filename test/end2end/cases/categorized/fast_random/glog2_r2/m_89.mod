var x{1..4} integer >= -2 <= 0;

maximize o: 2*x[1] + 3*x[2] + 5*x[3] + x[4];

s.t. c1: ((if x[1] != 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c2: ((x[1] > 0) ==> (x[4] < -1)) and ((x[4] < -1) ==> (x[3] >= 0));
s.t. c3: atleast 3 {i in 1..4} (x[i] >= -2);
s.t. c4: ((if x[1] != 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c5: atleast 3 {i in 1..4} (x[i] >= -2);
