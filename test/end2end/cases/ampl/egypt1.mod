
# =================================================
# EGYPT:  GAMS Egyptian Fertilizer Model
# =================================================

     # ORDERED PAIRS VERSION

     # Source:  "On the Development of a General Algebraic Modeling 
     # System in a Strategic Planning Environment" by Johannes Bisschop 
     # and Alexander Meeraus

# --------------------------------------
# Sets
# --------------------------------------

set center;                # Locations from which final product may be shipped

set port within center;    # Locations at which imports can be received

set plant within center;   # Locations of plants

set region;                # Demand regions

set unit;                  # Productive units

set proc;                  # Processes

set nutr;                  # Nutrients

set c_final;               # Final products (fertilizers)

set c_inter;               # Intermediate products

set c_ship within c_inter; # Intermediates for shipment

set c_raw;                 # Domestic raw materials and miscellaneous inputs

set commod := c_final union c_inter union c_raw;
			   # All commodities

# --------------------------------------
# Parameters
# --------------------------------------

param cf75 {region,c_final} >= 0; 
			   # Consumption of fertilizer 1974-75 (1000 tpy)

param fn {c_final,nutr} >= 0;
			   # Nutrient content of fertilizers

param cn75 {r in region, n in nutr} 
	      := sum {c in c_final} cf75[r,c] * fn[c,n];

			   # Consumption of nutrients 1974-75 (1000 tpy)

param road {region,center} >= 0;
			   # Road distances

param rail_half {plant,plant} >= 0;
param rail {p1 in plant, p2 in plant} := rail_half[p1,p2] + rail_half[p2,p1];

			   # Interplant rail distances (kms)

param impd_barg {plant} >= 0;
param impd_road {plant} >= 0;
			   # Import distances (kms) by barge and road

param tran_final {pl in plant, r in region} := 
	      if road[r,pl] > 0 then (.5 + .0144 * road[r,pl]) else 0;

param tran_import {r in region, po in port} :=
	      if road[r,po] > 0 then (.5 + .0144 * road[r,po]) else 0;

param tran_inter {p1 in plant, p2 in plant} :=
	      if rail[p1,p2] > 0 then (3.5 + .03 * rail[p1,p2]) else 0;

param tran_raw {pl in plant} :=
	    (if impd_barg[pl] > 0 then (1.0 + .0030 * impd_barg[pl]) else 0)
	  + (if impd_road[pl] > 0 then (0.5 + .0144 * impd_road[pl]) else 0);

			   # Transport cost (le per ton) for:
			   #   final products, imported final products,
			   #   interplant shipment, imported raw materials

param io {commod,proc};    # Input-output coefficients

param util {unit,proc} >= 0;
			   # Capacity utilization coefficients

param p_imp {commod} >= 0;
			   # Import Price (cif US$ per ton 1975)

param p_r {c_raw} >= 0;
param p_pr {plant,c_raw} >= 0;

param p_dom {pl in plant, c in c_raw} := 
	      if p_r[c] > 0 then p_r[c] else p_pr[pl,c];

			   # Domestic raw material prices

param dcap {plant,unit} >= 0;
			   # Design capacity of plants (t/day)

param icap {u in unit, pl in plant} := 0.33 * dcap[pl,u];
			   # Initial capacity of plants (t/day)

param exch := 0.4;         # Exchange rate

param util_pct := 0.85;    # Utilization percent for initial capacity

# --------------------------------------
# Derived sets of "possibilities"
# --------------------------------------

set m_pos := {u in unit, pl in plant: icap[u,pl] > 0};

			   # Set of (unit,plant) pairs for which there is
			   # initial capacity

set p_cap := {pr in proc, pl in plant: 
		      forall {u in unit: util[u,pr] > 0} (u,pl) in m_pos};

			   # Set of (proc,plant) pairs such that all units
			   # needed by the proc have some initial capacity
			   # at the plant

set p_except within proc cross plant;
			   # List of (pr,pl) such that process pr is 
			   # arbitrarily ruled out at plant pl

set p_pos := p_cap diff p_except;

			   # Set of possible (pr,pl) pairs
			   
set cp_pos := {c in commod, pl in plant: 
			       sum {pr in proc: (pr,pl) in p_pos} io[c,pr] > 0};

set cc_pos := {c in commod, pl in plant: 
			       sum {pr in proc: (pr,pl) in p_pos} io[c,pr] < 0};

set c_pos := cp_pos union cc_pos;

			   # Sets of commodity possibilities at each plant:
			   # production (cp_pos), consumption (cc_pos), 
			   # either production or consumption (c_pos)

# --------------------------------------
# Variables
# --------------------------------------

var Z {p_pos} >= 0;        # Z[pr,pl] is level of process pr at plant pl

var Xf {c in c_final, pl in plant, region: (c,pl) in cp_pos} >= 0;

			   # Xf[c,pl,r] is amount of final product c 
			   # shipped from plant pl to region r

var Xi {c in c_ship, p1 in plant, p2 in plant:
	      (c,p1) in cp_pos and (c,p2) in cc_pos} >= 0;

			   # Xi[c,p1,p2] is amount of intermediate c
			   # shipped from plant p1 to plant p2

var Vf {c_final,region,port} >= 0;
			   # Vf[c,r] is amount of final product c
			   # imported by region r from port po

var Vr {c in c_raw, pl in plant: (c,pl) in cc_pos} >= 0; 

			   # Vr[c,pl] is amount of raw material c
			   # imported for use at plant pl

var U {c in c_raw, pl in plant: (c,pl) in cc_pos} >= 0; 

			   # U[c,pl] is amount of raw material c
			   # purchased domestically for use at plant pl

var Psip;                  # Domestic recurrent cost
var Psil;                  # Transport cost
var Psii;                  # Import cost

# --------------------------------------
# Objective
# --------------------------------------

minimize Psi:  Psip + Psil + Psii;

# --------------------------------------
# Constraints
# --------------------------------------

subject to mbd {n in nutr, r in region}:

    sum {c in c_final} fn[c,n] * 
		(sum {po in port} Vf[c,r,po] +
		 sum {pl in plant: (c,pl) in cp_pos} Xf[c,pl,r])  >=  cn75[r,n];

			   # Total nutrients supplied to a region by all 
			   # final products (sum of imports plus internal 
			   # shipments from plants) must meet requirements
			   
subject to mbdb {c in c_final, r in region: cf75[r,c] > 0}:

    sum {po in port} Vf[c,r,po] +
    sum {pl in plant: (c,pl) in cp_pos} Xf[c,pl,r]  >=  cf75[r,c];

			   # Total of each final product supplied to each
			   # region (as in previous constraint) must meet
			   # requirements

subject to mb {c in commod, pl in plant}:

    sum {pr in proc: (pr,pl) in p_pos} io[c,pr] * Z[pr,pl]

   + ( if (c in c_ship) then sum {p2 in plant}
	  ( (if (c,p2) in cp_pos and (c,pl) in cc_pos then Xi[c,p2,pl])
	  - (if (c,p2) in cc_pos and (c,pl) in cp_pos then Xi[c,pl,p2]) ) )

   + ( if (c in c_raw and (c,pl) in cc_pos) then 
	 ( (if p_imp[c] > 0 then Vr[c,pl]) + (if p_dom[pl,c] > 0 then U[c,pl]) ) ) 

  >= if (c in c_final and (c,pl) in cp_pos) then sum {r in region} Xf[c,pl,r];

			   # For each commodity at each plant:  sum of
			   # (1) production or consumption at plant,
			   # (2) inter-plant shipments in or out,
			   # (3) import and domestic purchases (raw only)
			   # is >= 0 for raw materials and intermediates;
			   # is >= the total shipped for final products

subject to cc {(u,pl) in m_pos}:

    sum {pr in proc: (pr,pl) in p_pos} 
			   util[u,pr] * Z[pr,pl]  <=  util_pct * icap[u,pl];

			   # For each productive unit at each plant,
			   # the level of each process in that unit
			   # may not exceed the unit's capacity

subject to ap:

    Psip  =  sum {(c,pl) in cc_pos: c in c_raw} p_dom[pl,c] * U[c,pl];

			   # Psip is the cost of domestic raw materials,
			   # summed over all plants that consume them

subject to al:

    Psil  =  sum {c in c_final} (

		 sum {pl in plant, r in region: (c,pl) in cp_pos} 
					     tran_final[pl,r] * Xf[c,pl,r]

	       + sum {po in port, r in region} tran_import[r,po] * Vf[c,r,po] )

	   + sum {c in c_ship, p1 in plant, p2 in plant: 
			(c,p1) in cp_pos and (c,p2) in cc_pos} 
					   tran_inter[p1,p2] * Xi[c,p1,p2]

	   + sum {c in c_raw, pl in plant:
			(c,pl) in cc_pos and p_imp[c] > 0} 
						   tran_raw[pl] * Vr[c,pl];

			   # Total transport cost is sum of shipping costs for
			   # (1) all final products from all plants,
			   # (2) all imports of final products,
			   # (3) all intermediates shipped between plants,
			   # (4) all imports of raw materials

subject to ai:

    Psii / exch  =  sum {c in c_final, r in region, po in port} 
						       p_imp[c] * Vf[c,r,po]

		  + sum {c in c_raw, pl in plant: (c,pl) in cc_pos}
						       p_imp[c] * Vr[c,pl];

			   # Total import cost -- at exchange rate -- 
			   # is sum of import costs for final products 
			   # in each region and raw materials at each plant
