
# =================================================
# OIL:  GAMS Oil Refinery Model -- slice version
# =================================================

     # Source:  "Oil-Refinery Modeling with the
     # GAMS Language" by David Kendrick,
     # Alexander Meeraus and Jung Sun Suh

# --------------------------------------
# Sets
# --------------------------------------

set crude;                 # Types of crude oil

set m_raw;                 # Raw materials

set m_inter;               # Intermediates

set m_purch in m_inter;    # Purchased intermediates

set m_final;               # Final products

set proc;                  # Processes

set unit;                  # Productive units

set qual;                  # Quality attributes

set blend within {m_final, m_inter};

                           # (mf,mi) in this set implies that intermediate
                           # mi can be blended in making final product mf

# --------------------------------------
# Parameters
# --------------------------------------

param io {m_raw union m_inter, crude, proc} default 0;

                           # io[m,c,p] is input (if -) or output (if +)
                           # of material m, when using crude type c,
                           # per unit level of operation of process p

param util {unit,proc} logical default 0; 

                           # util[u,p] is 1 if process p uses unit u,
                           # 0 if process p does not use unit u

param cap {unit} >= 0 default 0;

                           # cap[u] is initial capacity of unit u

param price {crude union m_purch union m_final} > 0;

                           # price[m] is price per barrel of material m

param cost {proc} > 0;     # cost[p] is $ per barrel to operate process p

param attr_min {m_final,qual} >= 0 default 0;
param attr_max {m_final,qual} >= 0 default 0;

                           # lower and upper bounds on certain quality 
                           # attributes for the final products

param attr_both {m_inter,qual} >= 0 default 0;

param attr_crude {m_inter,crude,qual} >= 0 default 0;

     check {mi in m_inter, c in crude, q in qual}: 
               attr_crude[mi,c,q] = 0 or attr_both[mi,q] = 0;

param attr {mi in m_inter, c in crude, q in qual} 

     := if attr_crude[mi,c,q] > 0 then attr_crude[mi,c,q] else attr_both[mi,q];

                           # attr[mi,c,q] is amount of quality attribute q
                           # contributed to blend by intermediate mi that
                           # was derived from crude c

param purch_max {crude} > 0;

                           # upper bounds on purchases of crude oil

# --------------------------------------
# Variables
# --------------------------------------

var InCr {crude} >= 0;     # InCr[c] is purchase of crude oil c

var InInt {m_purch,crude} >= 0;  
                           # InInt[mi,c] is purchase of intermediate material
                           # mi derived from crude c


var LevPr {proc,crude} >= 0;
                           # LevPr[p,c] is level of process p using crude c

var LevBl {blend,crude} >= 0;
                           # LevBl[mf,mi,c] is level of intermediate mi,
                           # from crude c, blended into final product mf

var Out {m_final} >= 0;    # Out[mf] is amount of final product mf made

var Revenue;               # Revenue from all final product sales

var PurchCost;             # Cost of all purchased materials

var OperCost;              # Cost of all operations

# --------------------------------------
# Objective
# --------------------------------------

maximize net:  Revenue - PurchCost - OperCost;

                           # Maximize revenue from sales less all costs

# --------------------------------------
# Constraints
# --------------------------------------

subject to bal_crude {mr in m_raw, c in crude}:

  sum {p in proc} io[mr,c,p] * LevPr[p,c] + InCr[c] >= 0;

                           # Input of each crude to all processes (negative
                           # by convention) must be covered by purchases

subject to bal_inter {mi in m_inter, c in crude}:

  sum {p in proc} io[mi,c,p] * LevPr[p,c]
      + (if mi in m_purch then InInt[mi,c])

	   >= sum {(mf,mi) in blend} LevBl[mf,mi,c];

                           # Amount of each intermediate produced and 
                           # purchased must cover amount used in blending

subject to bal_final {mf in m_final}:

  Out[mf] = sum {(mf,mi) in blend, c in crude} LevBl[mf,mi,c];

                           # Volume of each final product must equal volume
                           # of all intermediates blended into it

subject to qual_min {mf in m_final, q in qual: attr_min[mf,q] <> 0}:

  sum {mi in m_inter, c in crude: (mf,mi) in blend}
        attr[mi,c,q] * LevBl[mf,mi,c] >= attr_min[mf,q] * Out[mf];

                           # Quality attributes of final products must
                           # not fall below specified minimums

subject to qual_max {mf in m_final, q in qual: attr_max[mf,q] <> 0}:

  sum {mi in m_inter, c in crude: (mf,mi) in blend}
        attr[mi,c,q] * LevBl[mf,mi,c] <= attr_max[mf,q] * Out[mf];

                           # Quality attributes of final products must
                           # not exceed specified maximums

subject to capacity {u in unit}:

  sum {p in proc} (util[u,p] * sum {c in crude} LevPr[p,c]) <= cap[u];

                           # Use of a productive unit by all processes
                           # must not exceed specified capacity

subject to purchase {c in crude}:  InCr[c] <= purch_max[c];

                           # Purchases of crude must not exceed specified
                           # maximum

subject to def_Revenue:

  Revenue = sum {mf in m_final} price[mf] * Out[mf];

                           # Total revenue is sum over final products of
                           # price times output

subject to def_PurchCost:

  PurchCost = sum {c in crude} price[c] * InCr[c] +
              sum {mp in m_purch, c in crude} price[mp] * InInt[mp,c];

                           # Total cost of purchases is cost of crudes plus
                           # cost of intermediates that are purchased

subject to def_OperCost:

  OperCost = sum {p in proc} (cost[p] * sum {c in crude} LevPr[p,c]);

                           # Total operating cost is sum over processes
                           # of cost times level for all crudes
