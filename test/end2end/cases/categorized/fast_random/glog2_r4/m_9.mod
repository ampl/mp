var x{1..3} integer >= 0 <= 2;

maximize o: 5*x[1] + x[2] - x[3];

s.t. c1: alldiff(x[3], x[1]);
s.t. c2: (x[2] < 2) ==> (x[1] <= 2) else (x[2] == 2);
s.t. c3: atmost 1 {i in 1..3} (x[i] >= 0);
s.t. c4: atmost 1 {i in 1..3} (x[i] >= 1);
s.t. c5: alldiff(x[3], x[1]);
s.t. c6: atmost 1 {i in 1..3} (x[i] >= 1);
