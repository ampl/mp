# ==>pies.mod

# A linear program with a variable rhs in the constraint system
# is expressed as a complementarity problem.
#
# LP:	min 	<c,x>
# 		Ax = q(p),
# 	s.t.	Bx = b,
# 		x >= 0
#
# where the prices p are the duals to the first constraint.
#
# MCP:	       A'p + B'v + c >= 0, x >= 0, comp.
# 	-Ax + q(p)           =  0, p free, comp.
# 	-Bx              + b =  0, v free, comp.
#
# Of course, the variables x and v are going to be
# split up further in the model.

# References:
#	William W. Hogan, "Energy Policy Models for Project Independence",
#		Computers & Operations Research (2), 1975.
#
#	N. Josephy, "A Newton Method for the PIES Energy Model",
#		Tech Report, Mathematics Research Center, UW-Madison, 1979.

set COMOD := {"C","L","H"};	# coal and light and heavy oil

set R;		# resources
set CREG;	# coal producing regions
set OREG;	# crude oil producing regions
set CTYP;	# increments of coal production
set OTYP;	# increments of oil production
set REFIN;	# refineries
set USERS;	# consumption regions

param rmax {R}; 	# maximum resource usage
param cmax {CREG,CTYP};	# coal prod. limits
param omax {OREG,OTYP};	# oil prod. limits
param rcost {REFIN};	# REFINing cost
param q0 {COMOD};	# base demand for commodities
param p0 {COMOD};	# base prices for commodities
param demand {COMOD,USERS};	# computed at optimality
param output {REFIN,COMOD};	# % output of light/heavy oil
param esub {COMOD,COMOD};	# cross-elasticities of substitution
param cruse {R,CREG,CTYP};	# resource use in coal prod
param oruse {R,OREG,OTYP};	# resource use in oil prod
param ccost {CREG,CTYP};	# coal prod. cost
param ocost {OREG,OTYP};	# oil prod. cost
param ctcost {CREG,USERS};
param otcost {OREG,REFIN};
param ltcost {REFIN,USERS};	# light oil trans costs
param htcost {REFIN,USERS};	# heavy oil trans costs

var C {CREG, CTYP};	# coal production
var O {OREG, OTYP};	# oil production
var Ct {CREG, USERS};	# coal transportation levels
var Ot {OREG, REFIN};	# crude oil transportation levels
var Lt {REFIN, USERS};	# light transportation levels
var Ht {REFIN, USERS};	# heavy transportation levels
var P {COMOD, USERS};	# commodity prices
var Mu {R} := 1;	# dual to ruse cons.; marginal utility
var Cv {CREG} := 1;	# dual to cmbal
var Ov {OREG} := 1;
var Lv {REFIN} := 1;
var Hv {REFIN} := 1;

s.t. delc {c in CREG, t in CTYP}:
  0 <= C[c,t] <= cmax[c,t] complements
  ccost[c,t] + (sum {res in R} cruse[res,c,t] * Mu[res]) - Cv[c];

s.t. delo {o in OREG, t in OTYP}:
  0 <= O[o,t] <= omax[o,t] complements
  ocost[o,t] + (sum {res in R} oruse[res,o,t] * Mu[res]) - Ov[o];

s.t. delct {c in CREG, u in USERS}:
  0 <= Ct[c,u] complements
  ctcost[c,u] + Cv[c] >= P["C",u];

s.t. delot {o in OREG, r in REFIN}:
  0 <= Ot[o,r] complements
  otcost[o,r] + rcost[r] + Ov[o] >=
     output[r,"L"] * Lv[r] + output[r,"H"] * Hv[r];

s.t. dellt {r in REFIN, u in USERS}:
  0 <= Lt[r,u] complements
  ltcost[r,u] + Lv[r] >= P["L",u];

s.t. delht {r in REFIN, u in USERS}:
  0 <= Ht[r,u] complements
  htcost[r,u] + Hv[r] >= P["H",u];

s.t. dembal {co in COMOD, u in USERS}:	# excess supply of product
  .1 <= P[co,u] complements
  (if co = "C" then sum {c in CREG} Ct[c,u]) +
  (if co = "L" then sum {r in REFIN} Lt[r,u]) +
  (if co = "H" then sum {r in REFIN} Ht[r,u]) >=
     q0[co] * prod {cc in COMOD} (P[cc,u]/p0[cc])**esub[co,cc];

s.t. cmbal {c in CREG}:			# coal material balance
  Cv[c] complements 
  sum {t in CTYP} C[c,t] = sum {u in USERS} Ct[c,u];
  
s.t. ombal {o in OREG}:			# oil material balance
  Ov[o] complements
  sum {t in OTYP} O[o,t] = sum {r in REFIN} Ot[o,r];

s.t. lmbal {r in REFIN}:		# light material balance
  Lv[r] complements
  sum {o in OREG} Ot[o,r] * output[r,"L"] = sum {u in USERS} Lt[r,u];

s.t. hmbal {r in REFIN}:		# heavy material balance
  Hv[r] complements
  sum {o in OREG} Ot[o,r] * output[r,"H"] = sum {u in USERS} Ht[r,u];

s.t. ruse {res in R}:			# resource use constraints
  0 <= Mu[res] complements
  rmax[res] >=
     sum {c in CREG, t in CTYP} C[c,t] * cruse[res,c,t] +
     sum {o in OREG, t in OTYP} O[o,t] * oruse[res,o,t];
