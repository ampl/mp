# From Bob Fourer's TOMS paper, June 1983


# A factory can manufacture some number of different products
# over the next T production periods.  Each product returns a
# characteristic estimated profit per unit, which varies from
# period to period.  The factory's size imposes a fixed upper
# limit on the total units manufactured per period.  Additionally,
# each product requires fixed characteristic amounts of certain
# raw materials per unit.
# 
# Limited quantities of raw materials must be stored now for use in 
# the next T periods.  Each raw material has a fixed characteristic
# storage cost per unit period.  Any material still unused after
# period T has a certain estimated remaining value.
# 
# What products should be manufactured in what periods to maximize
# total expected profit minus total storage costs, adjusted for
# the remaining value of any unused raw materials?



set P;			# Products
set R;			# Raw materials

param T > 0;		# number of production periods
param M > 0;		# Maximum production per period
param a{R,P} >= 0;	# units of raw material i to manufacture
			# 1 unit of product j
param b{R} >= 0;	# maximum initial stock of raw material i
param c{P,1..T};	# estimated profit per unit of product in period t
param d{R};		# storage cost per period per unit of raw material
param f{R};		# estimated remaining value per unit of raw material
			# after last period
var x{P,1..T} >= 0;	# units of product manufactured in period
var s{R,1..T+1} >= 0;	# stock of raw material at beginning of period

maximize profit:
	sum {t in 1..T} ( sum {j in P} c[j,t]*x[j,t] - sum {i in R} d[i]*s[i,t] )
	+ sum {i in R} f[i] * s[i,T+1];
		# total over all periods of estimated profit - storage costs
		# + value of raw materials left after last period

subject to
prd {t in 1..T}:
	sum {j in P} x[j,t] <= M;
		# production in each period less than maximum
stock1 {i in R}:
	s[i,1] <= b[i];
		# stock in period 1 less than maximum
stock {i in R, t in 1..T}:
	s[i,t+1] = s[i,t] - sum {j in P} a[i,j] * x[j,t];
		# stock in next period = present period - raw materials
