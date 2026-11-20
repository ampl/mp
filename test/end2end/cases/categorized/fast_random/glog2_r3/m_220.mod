var x{1..5} integer >= 0 <= 2;

maximize o: -x[1] - x[2] - 2*x[3] + 2*x[4] - 5*x[5];

s.t. c1: atmost 1 {i in 1..5} (x[i] >= 2);
s.t. c2: atleast 1 {i in 1..5} (x[i] >= 2);
s.t. c3: (if x[5] == 2 then x[2] else x[1]) >= 1;
