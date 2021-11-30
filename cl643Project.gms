$title eps-Constraint Method for Multiobjective Optimization (EPSCM,SEQ=319)

$onText
The eps-Constraint Method

This is a GAMS implementation of the augmented eps-constraint method
for generating the efficient (Pareto optimal, nondominated) solutions
in multiobjective problems. The eps-constraint method optimizes one of
the objective functions using the remaining objective functions as
constraints, varying their right hand side.

The generated optimal solutions proved to be efficient solutions of
the multiobjective problem under certain conditions.

The eps-constraint method consists of two phases:
1. Creation of the payoff table
2. Use the ranges from the payoff table in order to apply the method

The augmented eps-constraint uses lexicographic optimization in the
construction of the payoff table (in order to secure the Pareto
optimality of the individual optima) and a slightly modified objective
function in order to ensure the production of Pareto optimal (and not
weakly Pareto optimal) solutions. In addition, it performs early exit
from infeasible loops improving the performance of the algorithm in
multi-objective problems with several objective functions.

The algorithm can work also with MIP models. Actually the advantages
of the eps-constraint method over the weighting method are more
apparent for MIP problems where the non supported Pareto optimal
solutions can be produced.


$offText

$if not setEnv GMSPYTHONLIB $abort.noError Embedded code Python not ready to be used
$log --- Using Python library %sysEnv.GMSPYTHONLIB%

$inlineCom [ ]
$eolCom //

$stitle Example Model Definitions
Set
   i       'affected area' /A1*A11/
   j       'shelter' /S1*S20/
   k       'objective functions' /cost, time/;
$set min -1
$set max +1

Parameters
   dir(k) 'direction of the objective functions'
          /cost -1, time -1/
   d(i,j)  'Distance between affected area i and candidate shelter j'
   c(j)    'Capacity of the candidate shelter j'
    /S1 500
     S2 2000
     S3 500
     S4 2000
     S5 500
     S6 2000
     S7 2000
     S8 2000
     S9 500
     S10 500
     S11 2000
     S12 2000
     S13 2000
     S14 2000
     S15 2000
     S16 2000
     S17 500
     S18 2000
     S19 2000
     S20 500/
   fc 'Fuel cost';
   fc = 80;
   
Table d(i,j) 'matrix of distance between i and j'
    S1    S2    S3    S4    S5    S6    S7    S8    S9   S10   S11   S12   S13   S14   S15   S16   S17   S18   S19   S20 
A1  19    18    34    17    13    38    29    32    32    14    19    38    30    30    24    40    23    11    13    37
A2  38    28    10    25    18    37    36    25    14    15    23    11    26    15    13    32    19    12    28    26
A3  23    32    38    29    20    20    34    24    30    22    10    32    31    21    28    20    15    19    34    39
A4  15    16    32    31    31    31    27    11    26    35    40    18    30    24    17    28    15    26    23    29
A5  38    13    25    22    14    16    15    31    40    34    15    23    15    40    21    13    23    30    12    39
A6  40    19    27    21    32    10    17    11    30    11    13    26    13    14    28    38    12    22    18    17
A7  23    19    17    40    13    33    37    12    34    22    21    39    40    36    17    37    28    35    14    30
A8  13    23    24    11    30    25    10    26    24    26    16    22    15    29    19    35    24    32    18    18
A9  18    25    39    37    25    24    25    12    23    22    25    40    11    21    29    18    31    40    23    30
A10 22    12    26    38    34    38    15    35    35    30    20    19    27    15    18    28    31    26    26    31
A11 28    18    26    34    32    28    40    35    12    29    39    31    37    23    35    10    29    20    24    12;
 

Parameters  h(i)    'Number of victims in area i'
   /A1 1050 
    A2 740 
    A3 1100 
    A4 800 
    A5 650
    A6 550
    A7 450
    A8 450
    A9 900
    A10 426
    A11 490/;
      

Parameters

   m(i,j)  'Minimum acceptable distance between area i and shelter j';
   m(i,j) = 0.0011783*h(i);
  
   
Parameters
   Cv   'Capacity of vehicle'
   
   N       'Number of vehicles for evacuation'
   
   
   W       'Maximum allowed time for evacuating the victims from affected area i to shelter j'
   
   
   a       'Constant coefficient of transportation cost per kilometer per person per trip'
   
   
   b       'Wage per person per day for hiring staff to work in the shelter'
   
   
   Sv       'Ratio of the required staff per victim'
   
   
   T       'Duration of the disaster occurence'
   
   
   V       'Velocity of the vehicle using in evacuation process';
   
Cv = 8;
N = 50;
W = 72;
a = fc*0.025;
b = 400;
Sv = 0.2;
T = 6;
V = 32;

Variable
   O(k)  'Objective function variables'

Binary Variable
   X(i,j)    'Affected area i assigned to shelter j';

Positive Variable
   Y(i,j)    'Number of victims in area i that are assigned to shelter j';

Equation
   objcost   'objective for minimizing total cost'
   objtime    'objective for minimizing total time'


   Constraint1(i,j)   'Constraint3'
   Constraint2(i,j)   'Constraint4'
   Constraint3(j)   'Constraint5'
   Constraint4(i)   'Constraint6'
   Constraint5(i,j) 'Constraint7';

   
* Objective functions
objcost..   a*sum((i,j),d(i,j)*X(i,j)*h(i))+b*T*sum((i,j),Y(i,j)/Sv) =e= O('cost');

objtime..    sum((i,j),d(i,j)*X(i,j)*h(i))/(V*N*Cv) =e= O('time');


Constraint1(i,j)..  d(i,j)*X(i,j) =g= m(i,j);

Constraint2(i,j)..  d(i,j)*X(i,j)*h(i)/(V*N*Cv) =l= W;

Constraint3(j)..  sum(i,Y(i,j)) =l= c(j);

Constraint4(i)..  sum(j,Y(i,j)) =e= h(i);

Constraint5(i,j).. X(i,j) =g= 0;

Model example / all /;

$sTitle eps-constraint Method
Set
   k1(k)  'the first element of k'
   km1(k) 'all but the first elements of k'
   kk(k)  'active objective function in constraint allobj';

k1(k)$(ord(k) = 1) = yes;
km1(k)  = yes;
km1(k1) =  no;

Parameter
   rhs(k)    'right hand side of the constrained obj functions in eps-constraint'
   maxobj(k) 'maximum value from the payoff table'
   minobj(k) 'minimum value from the payoff table'
   numk(k)   'ordinal value of k starting with 1';

Scalar
   iter         'total number of iterations'
   infeas       'total number of infeasibilities'
   elapsed_time 'elapsed time for payoff and e-sonstraint'
   start        'start time'
   finish       'finish time';

Variable
   a_objval 'auxiliary variable for the objective function'
   obj      'auxiliary variable during the construction of the payoff table'
   sl(k)    'slack or surplus variables for the eps-constraints';

Positive Variable sl;

Equation
   con_obj(k) 'constrained objective functions'
   augm_obj   'augmented objective function to avoid weakly efficient solutions'
   allobj     'all the objective functions in one expression';

con_obj(km1).. O(km1) - dir(km1)*sl(km1) =e= rhs(km1);

* We optimize the first objective function and put the others as constraints
* the second term is for avoiding weakly efficient points

augm_obj..
   a_objval =e= sum(k1,dir(k1)*O(k1))
         + 1e-3*sum(km1,power(10,-(numk(km1) - 1))*sl(km1)/(maxobj(km1) - minobj(km1)));

allobj.. sum(kk, dir(kk)*O(kk)) =e= obj;

Model
   mod_payoff    / example, allobj            /
   mod_epsmethod / example, con_obj, augm_obj /;

Parameter payoff(k,k) 'payoff tables entries';

Alias (k,kp);

option optCr = 0, limRow = 0, limCol = 0, solPrint = off, solveLink = %solveLink.LoadLibrary%;

* Generate payoff table applying lexicographic optimization
loop(kp,
   kk(kp) = yes;
   repeat
      solve mod_payoff using mip maximizing obj;
      payoff(kp,kk) = O.l(kk);
      O.fx(kk) = O.l(kk); // freeze the value of the last objective optimized
      kk(k++1) = kk(k);   // cycle through the objective functions
   until kk(kp);
   kk(kp) = no;
*  release the fixed values of the objective functions for the new iteration
   O.up(k) =  inf;
   O.lo(k) = -inf;
   display objcost.l, objtime.l;
*   display Constraint1.l;
   display X.l;
   display Y.l;
   display payoff;
);
if(mod_payoff.modelStat <> %modelStat.optimal% and
   mod_payoff.modelStat <> %modelStat.integer Solution%,
   abort 'no optimal solution for mod_payoff';);

File fx / 2kp50_augmecon2_results.txt /;
put  fx ' PAYOFF TABLE'/;
loop(kp,
   loop(k, put payoff(kp,k):12:2;);
   put /;
);

minobj(k) = smin(kp,payoff(kp,k));
maxobj(k) = smax(kp,payoff(kp,k));
* gridpoints are calculated as the range (difference between max and min) of
* the 2nd objective function from the payoff table
$if not set gridpoints $set gridpoints 491
Set
   g         'grid points' / g0*g%gridpoints% /
   grid(k,g) 'grid';

Parameter
   gridrhs(k,g) 'RHS of eps-constraint at grid point'
   maxg(k)      'maximum point in grid for objective'
   posg(k)      'grid position of objective'
   firstOffMax  'some counters'
   lastZero     'some counters'
*  numk(k) 'ordinal value of k starting with 1'
   numg(g)      'ordinal value of g starting with 0'
   step(k)      'step of grid points in objective functions'
   jump(k)      'jumps in the grid points traversing';

lastZero = 1;
loop(km1,
   numk(km1) = lastZero;
   lastZero  = lastZero + 1;
);
numg(g) = ord(g) - 1;

grid(km1,g) = yes; // Here we could define different grid intervals for different objectives
maxg(km1)   = smax(grid(km1,g), numg(g));
step(km1)   = (maxobj(km1) - minobj(km1))/maxg(km1);
gridrhs(grid(km1,g))$(dir(km1) = -1) = maxobj(km1) - numg(g)/maxg(km1)*(maxobj(km1) - minobj(km1));
gridrhs(grid(km1,g))$(dir(km1) =  1) = minobj(km1) + numg(g)/maxg(km1)*(maxobj(km1) - minobj(km1));

put / ' Grid points' /;
loop(g,
   loop(km1, put gridrhs(km1,g):12:2;);
   put /;
);
put / 'Efficient solutions' /;

* Walk the grid points and take shortcuts if the model becomes infeasible or
* if the calculated slack variables are greater than the step size
posg(km1) = 0;
iter   = 0;
infeas = 0;
start  = jnow;

repeat
   rhs(km1) = sum(grid(km1,g)$(numg(g) = posg(km1)), gridrhs(km1,g));
   solve mod_epsmethod maximizing a_objval using mip;
   iter = iter + 1;
   if(mod_epsmethod.modelStat<>%modelStat.optimal% and
      mod_epsmethod.modelStat<>%modelStat.integer Solution%,
      infeas = infeas + 1; // not optimal is in this case infeasible
      put iter:5:0, '  infeasible' /;
      lastZero = 0;
      loop(km1$(posg(km1)  > 0 and lastZero = 0), lastZero = numk(km1));
      posg(km1)$(numk(km1) <= lastZero) = maxg(km1); // skip all solves for more demanding values of rhs(km1)
   else
      put iter:5:0;
      loop(k, put O.l(k):12:2;);
      jump(km1) = 1;
*     find the first off max (obj function that hasn't reach the final grid point).
*     If this obj.fun is k then assign jump for the 1..k-th objective functions
*     The jump is calculated for the innermost objective function (km=1)
      jump(km1)$(numk(km1) = 1) = 1 + floor(sl.L(km1)/step(km1));
      loop(km1$(jump(km1)  > 1), put '   jump';);
      put /;
   );
*  Proceed forward in the grid
   firstOffMax = 0;
   loop(km1$(posg(km1) < maxg(km1) and firstOffMax = 0),
      posg(km1)   = min((posg(km1) + jump(km1)),maxg(km1));
      firstOffMax = numk(km1);
   );
   posg(km1)$(numk(km1) < firstOffMax) = 0;
   abort$(iter > 1000) 'more than 1000 iterations, something seems to go wrong';
until sum(km1$(posg(km1) = maxg(km1)),1) = card(km1) and firstOffMax = 0;

finish = jnow;
elapsed_time = (finish - start)*60*60*24;

put /;
put 'Infeasibilities = ', infeas:5:0 /;
put 'Elapsed time: ',elapsed_time:10:2, ' seconds' /;