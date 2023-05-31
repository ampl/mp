# parameter values
param b default 1.25;
param p default 0.51;

# decision variables
var q1;
var q2 >= 0;
var w >= 0;

# objective
maximize ElogR:
    p*q1 + (1-p)*q2;

# conic constraints
s.t. T1: 1 + b*w >= exp( q1 );
s.t. T2: -1 + w +10*q2   <= -5 * q2 * exp( q1 / (q2*5) );
