## Economic order quantity
## MO-Book https://github.com/mobook/MO-book

param h default 0.75;      # cost of holding one item for one year 
param c default 500;       # cost of processing one order
param d default 10000;     # annual demand

# define variables for conic constraints
var x >= 0;
var y >= 0;

# conic constraint
s.t. q:
    x*y >= 1;

# linear objective
minimize eoq:
    h*x/2 + c*d*y;

