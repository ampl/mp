# Model DIST from June 1989 version of CSTR 133


#     This model determines a production and distribution plan to meet given
#demands for a set of goods.  The formulation is motivated by the experiences
#of a large producer in the United States.  The data are for a set of three
#products.


####  SHIPPING SETS AND PARAMETERS  ###

set whse 'warehouses';  # Locations from which demand is satisfied

set dctr 'distribution centers' within whse;
			# Locations from which product may be shipped

param sc 'shipping cost' {dctr,whse} >= 0;
			# Shipping costs, to whse from dctr, in $ / 100 lb

param huge 'largest shipping cost' > 0;
			# Largest cost allowed for a usable shipping route

param msr 'minimum size restriction' {dctr,whse} logical;
			# True indicates a minimum-size restriction on
			# direct shipments using this dctr --> whse route

param dsr 'direct shipment requirement' {dctr} >= 0;
			# Minimum total demand, in pallets, needed to
			# allow shipment on routes subject to the
			# minimum size restriction


###  PLANT SETS AND PARAMETERS  ###

set fact 'factories' within dctr;
			# Locations where product is manufactured

param rtmin 'regular-time total minimum' >= 0;
			# Lower limit on (average) total regular-time
			# crews employed at all factories

param rtmax 'regular-time total maximum' >= rtmin;
			# Upper limit on (average) total regular-time
			# crews employed at all factories

param otmin 'overtime total minimum' >= 0;
			# Lower limit on total overtime hours at all factories

param otmax 'overtime total maximum' >= otmin;
			# Upper limit on total overtime hours at all factories

param rmin 'regular-time minimums' {fact} >= 0;
			# Lower limits on (average) regular-time crews

param rmax 'regular-time maximums' {f in fact} >= rmin[f];
			# Upper limits on (average) regular-time crews

param omin 'overtime minimums' {fact} >= 0;
			# Lower limits on overtime hours

param omax 'overtime maximums' {f in fact} >= omin[f];
			# Upper limits on overtime hours

param hd 'hours per day' {fact} >= 0;
			# Regular-time hours per working day

param dp 'days in period' {fact} > 0;
			# Working days in the current planning period


###  PRODUCT SETS AND PARAMETERS  ###

set prod 'products';     # Elements of the product group

param wt 'weight' {prod} > 0;
			# Weight in 100 lb / 1000 cases

param cpp 'cases per pallet' {prod} > 0;
			# Cases of product per shipping pallet

param tc 'transshipment cost' {prod} >= 0;
			# Transshipment cost in $ / 1000 cases

param pt 'production time' {prod,fact} >= 0;
			# Crew-hours to produce 1000 cases

param rpc 'regular-time production cost' {prod,fact} >= 0;
			# Cost of production on regular time,
			# in $ / 1000 cases

param opc 'overtime production cost' {prod,fact} >= 0;
			# Cost of production on overtime, in $ / 1000 cases


###  DEMAND SETS AND PARAMETERS  ###

param dt 'total demand' {prod} >= 0;
			# Total demands for products, in 1000s

param ds 'demand shares' {prod,whse} >= 0.0, <= 1.0;

			# Historical demand data, from which each
			# warehouse's share of total demand is deduced

param dstot {p in prod} := sum {w in whse} ds[p,w];

			# Total of demand shares; should be 1, but often isn't

param dem 'demand' {p in prod, w in whse} := dt[p] * ds[p,w] / dstot[p];

			# Projected demands to be satisfied, in 1000s

set rt 'shipping routes available' :=

 {d in dctr, w in whse:
	 d <> w  and  sc[d,w] < huge  and
	 (w in dctr or sum {p in prod} dem[p,w] > 0)  and
	 not (msr[d,w] and sum {p in prod} 1000*dem[p,w]/cpp[p] < dsr[d]) };

			# List of ordered pairs that represent routes
			# on which shipments are allowed


###  VARIABLES  ###

var Rprd 'regular-time production' {prod,fact} >= 0;
			# Regular-time production of each product
			# at each factory, in 1000s of cases

var Oprd 'overtime production' {prod,fact} >= 0;
			# Overtime production of each product
			# at each factory, in 1000s of cases

var Ship 'shipments' {prod,rt} >= 0;
			# Shipments of each product on each allowed route,
			# in 1000s of cases

var Trans 'transshipments' {prod,dctr} >= 0;
			# Transshipments of each product at each
			# distribution center, in 1000s of cases



###  OBJECTIVE  ###

minimize cost:  sum {p in prod, f in fact} rpc[p,f] * Rprd[p,f] +
		sum {p in prod, f in fact} opc[p,f] * Oprd[p,f] +
		sum {p in prod, (d,w) in rt} sc[d,w] * wt[p] * Ship[p,d,w] +
		sum {p in prod, d in dctr} tc[p] * Trans[p,d];

			# Total cost:  regular production, overtime
			# production, shipping, and transshipment


###  CONSTRAINTS  ###

rtlim 'regular-time total limits':

    rtmin <= sum {p in prod, f in fact}
			(pt[p,f] * Rprd[p,f]) / (dp[f] * hd[f]) <= rtmax;

			# Total crews must lie between limits

otlim 'overtime total limits':

    otmin <= sum {p in prod, f in fact} pt[p,f] * Oprd[p,f] <= otmax;

			# Total overtime must lie between limits

rlim 'regular-time limits' {f in fact}:

    rmin[f] <= sum {p in prod}
			(pt[p,f] * Rprd[p,f]) / (dp[f] * hd[f]) <= rmax[f];

			# Crews at each factory must lie between limits

olim 'overtime limits' {f in fact}:

    omin[f] <= sum {p in prod} pt[p,f] * Oprd[p,f] <= omax[f];

			# Overtime at each factory must lie between limits

noRprd 'no regular production' {p in prod, f in fact: rpc[p,f] = 0}:

    Rprd[p,f] = 0;

noOprd 'no overtime production' {p in prod, f in fact: opc[p,f] = 0}:

    Oprd[p,f] = 0;      # Do not produce where specified cost is zero

bal 'material balance' {p in prod, w in whse}:

    sum {(v,w) in rt}
       Ship [p,v,w] + (if w in fact then Rprd[p,w] + Oprd[p,w]) =

    dem[p,w] + (if w in dctr then sum {(w,v) in rt} Ship[p,w,v]);

			# Demand is satisfied by shipment into warehouse
			# plus production (if it is a factory)
			# minus shipment out (if it is a distn. center)

trdef 'transshipment definition' {p in prod, d in dctr}:

    Trans[p,d] >= sum {(d,w) in rt} Ship [p,d,w] -
		  (if d in fact then Rprd[p,d] + Oprd[p,d]);

			# Transshipment at a distribution center is
			# shipments out less production (if any)
