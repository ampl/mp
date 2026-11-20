var x{1..3} integer >= -2 <= 1;
minimize o: (51)*x[1] + (-333)*x[2] + (2000)*x[3];
s.t. c1: atleast 3 {i in 1..3} (x[i] >= 1);
s.t. c2: ((x[3] <= -2) ==> (x[2] <= 1)) and ((x[2] <= 1) ==> (x[1] >= -2));
s.t. c3: (if x[2] < -2 then x[3] else x[3]) == -2;
