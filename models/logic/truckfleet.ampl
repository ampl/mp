# Truck configuration scheduling model.
# Adapted from truckfleet.cpp

include cp.ampl;
suffix priority IN;

param NumTruckConfigs integer > 0;
set TruckConfigs = 0..NumTruckConfigs - 1;

param NumOrders integer > 0;
set Orders = 0..NumOrders - 1;

param NumCustomers integer > 0;
set Customers = 0..NumCustomers - 1;

param NumTrips integer > 0;
set Trips = 0..NumTrips - 1;

param NumProductTypes integer > 0;
param MaxReconfigCost >= 0;

param MaxTruckConfigLoad{TruckConfigs} integer >= 0;
param TruckCost{TruckConfigs} integer >= 0;
param CustomerOfOrder{Orders} integer >= 0 <= NumCustomers - 1;
param Volumes{Orders} integer >= 0;
param ProductTypes{Orders} integer >= 0 <= NumProductTypes - 1;

# The cost for configuring the truck from one configuration to another
# depending on both configurations.
set ReconfigCosts dimen 3;

set AllowedContainerConfigs{0..NumProductTypes - 1};

var truckConfigs{Trips} integer >= 0 <= NumTruckConfigs - 1;
var orderTrip{Orders} integer >= 0 <= NumTrips - 1 suffix priority 1;
var load{Trips} integer >= 0 <= max{c in TruckConfigs} MaxTruckConfigLoad[c];
var customerOfTruck{Trips} integer >= 0 <= NumCustomers - 1;
var configOfContainer{Orders} integer >= 0 <= NumTruckConfigs - 1;
var reconfigCost{0..NumTrips - 2} integer >= 0 <= MaxReconfigCost;

minimize totalCost:
  sum{t in Trips}
    (if load[t] != 0 then element({c in TruckConfigs} TruckCost[c], truckConfigs[t])) +
  sum{t in 0..NumTrips - 2} reconfigCost[t];

minimize tripCount: count{t in 0..NumTrips - 1} (load[t] > 0);

s.t. reconfig{t in 0..NumTrips - 2}:
  in_relation(truckConfigs[t], truckConfigs[t + 1], reconfigCost[t],
              {(src, dst, cost) in ReconfigCosts} (src, dst, cost));

s.t. loadDef{t in Trips}:
  load[t] = sum{o in Orders} if orderTrip[o] = t then Volumes[o];

# Truck load cannot exceed max. load for truck's config.
s.t. satisfyMaxLoad{t in Trips}:
  load[t] <= element({c in TruckConfigs} MaxTruckConfigLoad[c], truckConfigs[t]);

s.t. configOfContainerDef{o in Orders}:
  configOfContainer[o] = element({t in Trips} truckConfigs[t], orderTrip[o]);

s.t. restrictConfigOfContainer{o in Orders}:
  in_relation(configOfContainer[o], {c in AllowedContainerConfigs[ProductTypes[o]]} c);

s.t. oneCustomerPerTruck{o in Orders}:
  element({t in Trips} customerOfTruck[t], orderTrip[o]) = CustomerOfOrder[o];

s.t. unusedTracksAtEnd{t in 1..NumTrips - 1}:
  load[t - 1] > 0 or load[t] = 0;

# Dominance: the non used truck keep the last used configuration.
s.t. loadFirstTruck: (load[0] > 0);
s.t. unusedKeepsLastConfig{t in 1..NumTrips - 1}:
  load[t] > 0 or truckConfigs[t] == truckConfigs[t - 1];

# Dominance: regroup deliveries with same configuration.
s.t. regroup{t in NumTrips - 2 .. 1}:
  truckConfigs[t] = truckConfigs[t - 1] or
  forall{u in t + 1 .. NumTrips - 1} truckConfigs[u] != truckConfigs[t - 1];

data;

param NumTruckConfigs := 7;
param NumOrders := 21;
param NumCustomers := 3;
param NumTrips := 15;
param NumProductTypes := 3;
param MaxReconfigCost := 1000;

param:
  MaxTruckConfigLoad TruckCost :=
0         11             2
1         11             2
2         11             2
3         11             3
4         10             3
5         10             3
6         10             4;

param:
   CustomerOfOrder Volumes ProductTypes :=
 0        0           3          1
 1        0           4          2
 2        0           3          0
 3        0           2          1
 4        0           5          1
 5        0           4          1
 6        0          11          0
 7        1           4          0
 8        1           5          0
 9        1           2          0
10        1           4          2
11        1           7          2
12        1           3          2
13        1           5          0
14        1           2          2
15        2           5          1
16        2           6          0
17        2          11          2
18        2           1          0
19        2           6          0
20        2           3          0;

set ReconfigCosts :=
  0 0 0
  0 1 0
  0 2 0
  0 3 10
  0 4 10
  0 5 10
  0 6 15
  1 0 0
  1 1 0
  1 2 0
  1 3 10
  1 4 10
  1 5 10
  1 6 15
  2 0 0
  2 1 0
  2 2 0
  2 3 10
  2 4 10
  2 5 10
  2 6 15
  3 0 3
  3 1 3
  3 2 3
  3 3 0
  3 4 10
  3 5 10
  3 6 15
  4 0 3
  4 1 3
  4 2 3
  4 3 10
  4 4 0
  4 5 10
  4 6 15
  5 0 3
  5 1 3
  5 2 3
  5 3 10
  5 4 10
  5 5 0
  5 6 15
  6 0 3
  6 1 3
  6 2 3
  6 3 10
  6 4 10
  6 5 10
  6 6 0;

set AllowedContainerConfigs[0] := 0 3 4 6;
set AllowedContainerConfigs[1] := 1 3 5 6;
set AllowedContainerConfigs[2] := 2 4 5 6;

option solver 'ilogcp';
option ilogcp_options 'outlev=1 timelimit=20 logperiod=50000';
solve;

print 'Configuration cost:', totalCost, 'Number of Trips:', tripCount;
print{t in Trips: load[t] > 0}: 'Trip', t & ': Config =', truckConfigs[t],
  'Items =', {o in Orders: orderTrip[o] = t} ('<' & o & ',' & ProductTypes[o] & ',' & Volumes[o] & '>');
