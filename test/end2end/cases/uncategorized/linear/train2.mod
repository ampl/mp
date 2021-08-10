# Model TRAIN from June 1989 version of CSTR 133


#     Given a day's schedule, this model allocates passenger cars to
#trains so as to minimize either the number of cars required or the
#number of car-miles run (Fourer, Gertler and Simkowitz 1977, 1978).  The
#data represent a hypothetical schedule and demands for service between
#Washington, Philadelphia, New York and Boston.


####  SCHEDULE SETS AND PARAMETERS  ###

set cities;

set links within {c1 in cities, c2 in cities: c1 <> c2};

                        # Set of cities, and set of intercity links

param last > 0 integer; # Number of time intervals in a day

set times := 1..last;   # Set of time intervals in a day

set schedule within
      {c1 in cities, t1 in times,
       c2 in cities, t2 in times: (c1,c2) in links};

                        # Member (c1,t1,c2,t2) of this set represents
                        # a train that leaves city c1 at time t1
                        # and arrives in city c2 at time t2


###  DEMAND PARAMETERS  ###

param section > 0 integer;
                        # Maximum number of cars in one section of a train

param demand {schedule} > 0;
                        # For each scheduled train:
                        # the smallest number of cars that
                        # can meet demand for the train

param low {(c1,t1,c2,t2) in schedule} := ceil(demand[c1,t1,c2,t2]);

                        # Minimum number of cars needed to meet demand

param high {(c1,t1,c2,t2) in schedule}

   := max (2, min (ceil(2*demand[c1,t1,c2,t2]),
                   section*ceil(demand[c1,t1,c2,t2]/section) ));

                        # Maximum number of cars allowed on a train:
                        # 2 if demand is for less than one car;
                        # otherwise, lesser of
                        # number of cars needed to hold twice the demand, and
                        # number of cars in minimum number of sections needed


###  DISTANCE PARAMETERS  ###

param dist_table {links} >= 0 default 0.0;

param distance {(c1,c2) in links} > 0
   := if dist_table[c1,c2] > 0 then dist_table[c1,c2] else dist_table[c2,c1];

                        # Inter-city distances: distance[c1,c2] is miles
                        # between city c1 and city c2


###  VARIABLES  ###

var U 'cars stored' {cities,times} >= 0;
                        # u[c,t] is the number of unused cars stored
                        # at city c in the interval beginning at time t

var X 'cars in train' {schedule} >= 0;
                        # x[c1,t1,c2,t2] is the number of cars assigned to
                        # the scheduled train that leaves c1 at t1 and
                        # arrives in c2 at t2



###  OBJECTIVES  ###

minimize cars:
       sum {c in cities} U[c,last] +
       sum {(c1,t1,c2,t2) in schedule: t2 < t1} X[c1,t1,c2,t2];

                        # Number of cars in the system:
                        # sum of unused cars and cars in trains during
                        # the last time interval of the day

minimize miles:
       sum {(c1,t1,c2,t2) in schedule} distance[c1,c2] * X[c1,t1,c2,t2];

                        # Total car-miles run by all scheduled trains in a day


###  CONSTRAINTS  ###

account {c in cities, t in times}:

  U[c,t] = U[c, if t > 1 then t-1 else last] +

      sum {(c1,t1,c,t) in schedule} X[c1,t1,c,t] -
      sum {(c,t,c2,t2) in schedule} X[c,t,c2,t2];

                        # For every city and time:
                        # unused cars in the present interval must equal
                        # unused cars in the previous interval,
                        # plus cars just arriving in trains,
                        # minus cars just leaving in trains

satisfy {(c1,t1,c2,t2) in schedule}:

       low[c1,t1,c2,t2] <= X[c1,t1,c2,t2] <= high[c1,t1,c2,t2];

                        # For each scheduled train:
                        # number of cars must meet demand,
                        # but must not be so great that unnecessary
                        # sections are run

