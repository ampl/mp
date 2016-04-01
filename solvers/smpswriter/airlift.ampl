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

set AircraftType;
set Routes;

# The maximum number of flight hours for aircraft of type i available during
# the month.
param MaxFlightHours{i in AircraftType} >= 0;

# The number of flight hours required for an aircraft of type i to complete
# one flight of route j.
param FlightHours{i in AircraftType, j in Routes} >= 0;

# The number of flight hours required for aircraft of type i to fly
# route k, after having been switched from route j. Note that an increase of
# flights_increase[i, j, k] flights for route k results in the cancellation of
# (SwitchedFlightHours[i, j, k] / FlightHours[i, j]) * flights_increase[i, j, k]
# flights for route j, since "k flights" and "j flights" are not necessarily
# equal units.
param SwitchedFlightHours{i in AircraftType, j in Routes, k in Routes} >= 0;

# The carrying capacity (in tons) of a single flight of an aircraft of type i,
# flying route j.
param Capacity{i in AircraftType, j in Routes} >= 0;

# The number of flights originally planned for route j using aircraft of type i.
var flights{i in AircraftType, j in Routes} >= 0;

# The first-stage constraint.
s.t. flight_hours{i in AircraftType}:
  sum{j in Routes} FlightHours[i, j] * flights[i, j] <= MaxFlightHours[i];

# Increase in the number of flights for route k flown by aircraft type i,
# because of being switched from route j.
var flights_increase{i in AircraftType, j in Routes, k in Routes} >= 0;

# The second-stage constraint: we cannot switch away more flight hours from
# aircraft of type i and from route j, than we have originally scheduled for
# such.
s.t. switch_flight_hours{i in AircraftType, j in Routes}:
  sum{k in Routes: k != j}
    SwitchedFlightHours[i, j, j] * flights_increase[i, j, k]
      <= FlightHours[i, j] * flights[i, j]

# The recourse constraint that the demand for each route must be met.
# TODO
