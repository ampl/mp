###  SHIPPING SETS AND PARAMETERS  ###
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
			# Lower limit on total overtime hours 
			# at all factories
param otmax 'overtime total maximum' >= otmin;
			# Upper limit on total overtime hours 
			# at all factories
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

set prod 'products';    # Elements of the product group
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
param ds 'demand shares' {p in prod, w in whse} >= 0, <= 1;
			# Historical demand data, from which each
			# warehouse's share of total demand is deduced
param dstot {p in prod} := sum {w in whse} ds[p,w];
			# Total of demand shares;
			# should be 1, but often isn't
param dem 'demand' {p in prod, w in whse} := dt[p] * ds[p,w] / dstot[p];
			# Projected demands to be satisfied, in 1000s
set rt 'shipping routes available' :=
 {d in dctr, w in whse: 
	 d <> w  and  sc[d,w] < huge  and
	 (w in dctr or sum {p in prod} dem[p,w] > 0)  and
	 not (msr[d,w] and sum {p in prod} 1000*dem[p,w]/cpp[p] < dsr[d]) };
			# List of ordered pairs that represent routes
			# on which shipments are allowed

###  OBJECTIVE  ###

minimize cost: to_come; # Total cost: regular production, overtime
			# production, shipping, and transshipment

###  NODES (CONSTRAINTS)  ###

node RT:  rtmin <= net_out <= rtmax;
			# Source of all regular-time crews allocated
node OT:  otmin <= net_out <= otmax;
			# Source of all overtime hours allocated
node P_RT {fact};       # Sources of regular-time crews at factories
node P_OT {fact};       # Sources of overtime hours at factories
node M {prod,fact};     # Sources of manufacturing:
			# send to factory's W node for local demand;
			# send to factory's D node for distribution
node D {prod,dctr};     # Sources of distribution:
			# receive transshipped goods from center's W node;
			# receive manufactured goods from center's M node;
			# send to W nodes elsewhere
node W {p in prod, w in whse}:  net_in = dem[p,w];
			# Locations of warehousing:
			# receive from D nodes and local M node (if any),
			# to satisfy local demand;
			# send to local D node (if any) for transshipment

###  ARCS (VARIABLES)  ###

arc Work_RT {f in fact}
      from RT  to P_RT[f]  >= rmin[f],  <= rmax[f];  
			# Regular-time crews allocated to each factory
arc Work_OT {f in fact}
      from OT  to P_OT[f]  >= omin[f],  <= omax[f];  
			# Overtime hours allocated to each factory
arc Manu_RT {p in prod, f in fact: rpc[p,f] <> 0} >= 0
      from P_RT[f]  to M[p,f] (dp[f] * hd[f] / pt[p,f])
      obj cost (rpc[p,f] * dp[f] * hd[f] / pt[p,f]);
			# Regular-time crews allocated to 
			# manufacture of each product at each factory
arc Manu_OT {p in prod, f in fact: opc[p,f] <> 0} >= 0
      from P_OT[f]  to M[p,f] (1 / pt[p,f])  obj cost (opc[p,f] / pt[p,f]);
			# Overtime hours allocated to 
			# manufacture of each product at each factory
arc Prod_L {p in prod, f in fact} >= 0
      from M[p,f]  to W[p,f];
			# Manufacture of each product at each factory
			# to satisfy local demand, in 1000s of units
arc Prod_D {p in prod, f in fact} >= 0
      from M[p,f]  to D[p,f];
			# Manufacture of each product at each factory,
			# for distribution elsewhere, in 1000s of units
arc Ship {p in prod, (d,w) in rt} >= 0
      from D[p,d]  to W[p,w] obj cost (sc[d,w] * wt[p]);
			# Shipments of each product on each allowed route
arc Trans {p in prod, d in dctr} >= 0
      from W[p,d]  to D[p,d]  obj cost (tc[p]);
			# Transshipments of each product at each
			# distribution center
