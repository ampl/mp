
#  ==============================================
#  TRAIN:  Allocating railroad passenger cars
#  ==============================================

#  ----------------------------------------------
#  Data
#  ----------------------------------------------

set cities;             # Set of cities

param last > 0 integer; # Number of intervals into which
                        # a schedule-period is divided

set times := 1..last;   # Set of intervals into which
                        # a schedule-period is divided

set schedule within
	{c1 in cities, t1 in times, c2 in cities, t2 in times: c1 != c2};

			# Member (c1,t1,c2,t2) of this set represents
			# a train that leaves city c1 at time t1 
			# and arrives at city c2 at time t2

param section > 0 integer;
                        # Maximum number of cars in one section of a train

param demand {schedule} > 0 ; # was integer;
                        # For each scheduled train:
                        # the smallest number of cars that
                        # can meet demand for the train

param distance {(c1,c2) in setof {(a,b,c,d) in schedule}(a,c)} > 0;
			# Inter-city distances: distance[c1,c2] is miles
			# between city c1 and city c2

#  ----------------------------------------------
#  Variables
#  ----------------------------------------------

var U 'cars stored' {cities,times} >= 0;
			# u[c,t] is number of unused cars stored 
			# at city c in the interval beginning at time t

var X 'cars in train' {schedule} >= 0;
			# x[c1,t1,c2,t2] is number of cars assigned to
			# scheduled train that leaves c1 at t1 and
			# arrives c2 at t2

#  ----------------------------------------------
#  Objectives
#  ----------------------------------------------

minimize cars: 
       sum {c in cities} U[c,last] +
       sum {(c1,t1,c2,t2) in schedule: t2 < t1} X[c1,t1,c2,t2];

			# Number of cars in the system: 
			# sum of unused cars and cars in trains during
			# the last interval of the schedule-period

minimize miles: 
       sum {(c1,t1,c2,t2) in schedule} distance[c1,c2] * X[c1,t1,c2,t2];

			# Total car-miles run by all scheduled trains
			# in one schedule-period

#  ----------------------------------------------
#  Constraints
#  ----------------------------------------------

account {c in cities, t in times}:

  U[c,t] = U[c,if t > 1 then t-1 else last] +
           sum {(c1,t1,c,t) in schedule} X[c1,t1,c,t] -
           sum {(c,t,c2,t2) in schedule} X[c,t,c2,t2];

			# For every city and time: 
			# unused cars in present interval must equal
			# unused cars in previous interval, 
			# plus cars just arriving in trains,
			# minus cars just leaving in trains

satisfy {(c1,t1,c2,t2) in schedule}:

       demand[c1,t1,c2,t2] <= X[c1,t1,c2,t2] 
			   <= section * ceil(demand[c1,t1,c2,t2]/section);

			# For each scheduled train:
			# number of cars must meet demand, 
			# but must not be so great that unnecessary
			# sections are run
