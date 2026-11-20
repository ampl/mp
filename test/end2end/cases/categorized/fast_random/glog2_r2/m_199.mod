var x{1..5} integer >= 0 <= 4;

minimize o: x[1] + 4*x[2] + 3*x[3] + x[4] + 2*x[5];

s.t. c1: atleast 3 {i in 1..5} (x[i] >= 4);
s.t. c2: (x[5] != 2) ==> (x[5] == 4);
s.t. c3: exists {i in 1..5} x[i] == 0;
s.t. c4: ((if x[4] > 2 then 1 else 0) + (if x[5] >= 2 then 1 else 0)) >= 1;
s.t. c5: atleast 3 {i in 1..5} (x[i] >= 4);
s.t. c6: ((if x[4] > 2 then 1 else 0) + (if x[5] >= 2 then 1 else 0)) >= 1;
