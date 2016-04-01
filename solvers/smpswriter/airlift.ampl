# Airlift operation scheduling.
#
# References:
#
# * K. A. Ariyawansa and A. J. Felt. On a new collection of stochastic linear
#   programming test problems. Technical Report 4, Department of Mathematics,
#   Washington State University, Pullman, WA 99164, Apr. 2001.
#
# * J. L. Midler and R. D. Wollmer. Stochastic programming models for
#   scheduling airlift operations. Naval Research Logistics Quarterly,
#   16:315-330, 1969.
#
# AMPL coding by Victor Zverovich.

set AircraftTypes;
set Routes;

# The maximum number of flight hours for aircraft of type i available during
# the month (F_i).
param MaxHours{i in AircraftTypes} >= 0;

# The number of flight hours required for an aircraft of type i to complete
# one flight of route j (a_{ij}).
param Hours{i in AircraftTypes, j in Routes} >= 0;

# The number of flight hours required for aircraft of type i to fly
# route k, after having been switched from route j. Note that an increase of
# switched_flights[i, j, k] flights for route k results in the cancellation of
# (SwitchedFlightHours[i, j, k] / Hours[i, j]) * switched_flights[i, j, k]
# flights for route j, since "k flights" and "j flights" are not necessarily
# equal units (a_{ijk}).
param SwitchedHours{i in AircraftTypes, j in Routes, k in Routes} >= 0;

# The carrying capacity (in tons) of a single flight of an aircraft of type i,
# flying route j (b_{ij}).
param Capacity{i in AircraftTypes, j in Routes} >= 0;

# The cost for aircraft type i to be initially assigned and fly one flight of
# route j (c_{ij}).
param AssignCost{i in AircraftTypes, j in Routes} >= 0;

# The cost for aircraft type i to fly one flight of route k, after having been
# initially assigned route j (c_{ijk}).
param SwitchCost{i in AircraftTypes, j in Routes, k in Routes} >= 0;

# The cost per ton of commercially contracted transport on route j (c^+_j).
param ContractedCost{j in Routes} >= 0;

# The cost per ton of unused capacity on route j (c^-_j).
param UnusedCost{j in Routes} >= 0;

# The set of scenarios.
set Scen;

# The probability of scenario s.
param P{s in Scen} >= 0 <= 1;

# The random variable (parameter) representing the demand for route j (d_j).
param Demand{j in Routes} >= 0;

# The number of flights originally planned for route j using aircraft of
# type i (x_{ij}).
var flights{i in AircraftTypes, j in Routes} >= 0;

# Increase in the number of flights for route k flown by aircraft type i,
# because of being switched from route j (x_{ijk}).
var switched_flights{i in AircraftTypes, j in Routes, k in Routes} >= 0;

# The load originally scheduled to be carried on route j (i.e. the "best guess"
# of the demand).
var load{j in Routes} = sum{i in AircraftTypes} Capacity[i, j] * flights[i, j];

# The total carrying capacity switched away from route j in the recourse action.
var capacity_out{j in Routes} =
  sum{i in AircraftTypes, k in Routes: k != j}
    Capacity[i, j] * (SwitchedHours[i, j, k] / Hours[i, j]) *
      switched_flights[i, j, k];

# The carrying capacity switched to route j.
var capacity_in{j in Routes} =
  sum{i in AircraftTypes, k in Routes: k != j}
    Capacity[i, j] * switched_flights[i, k, j];

# The demand for route j which is contracted commercially in the
# recourse (y^+_j).
var contracted{j in Routes} >= 0;

# The unused capacity assigned to route j (y^-_j),
var unused{j in Routes} >= 0;

minimize expected_cost:
  sum{i in AircraftTypes, j in Routes} AssignCost[i, j] * flights[i, j] +
  sum{s in Scen} P[s] * (
    sum{i in AircraftTypes, j in Routes, k in Routes: k != j}
      (SwitchCost[i, j, k] -
        AssignCost[i, j] * (SwitchedHours[i, j, k] / Hours[i, j])) *
      switched_flights[i, j, k] +
    sum{j in Routes} ContractedCost[j] * contracted[j] +
    sum{j in Routes} UnusedCost[j] * unused[j]);

# The first-stage constraint.
s.t. flight_hours{i in AircraftTypes}:
  sum{j in Routes} Hours[i, j] * flights[i, j] <= MaxHours[i];

# The second-stage constraint: we cannot switch away more flight hours from
# aircraft of type i and from route j, than we have originally scheduled for
# such.
s.t. switch_flight_hours{i in AircraftTypes, j in Routes}:
  sum{k in Routes: k != j} SwitchedHours[i, j, j] * switched_flights[i, j, k]
      <= Hours[i, j] * flights[i, j];

# The recourse constraint that the demand for each route must be met.
s.t. satisfy_demand{j in Routes}:
  load[j] - capacity_out[j] + capacity_in[j] + contracted[j] - unused[j]
    = Demand[j];
