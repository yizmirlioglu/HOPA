# HOPA
# Hybrid Oriented Point Algebra

This repository includes the following items:

ASP folder:  ASP programs for HOPA

Scenario:  ASP encoding of the scenarios

HOPA_Instances: Consistent and inconsistent HOPA benchmark problem instances used in experiments

OPRA_Instances: Consistent and inconsistent OPRA benchmark problem instances used in experiments

Scripts: Python scripts to generate random benchmark problem instances and Python scripts to automatically test these instances

Technical Appendix.pdf : Proof of the Theorems, Full ASP encoding of HOPA

and also Clingo executable files to run the ASP solver and get the output.


Scenarios and experiments were performed on Ubuntu with Clingo 5.4.0  

For Windows OS, you should use the Clingo executable files in the "clingo-5.4.0-Windows64.zip" file. For this, you should replace executable files in the main folder with the ones in "clingo-5.4.0-Windows64.zip"


ASP PROGRAMS:

hopa_main.lp :  The main ASP program for checking consistency of HOPA constraints

opra_consistency.lp :  The simplified version of the main ASP program which is used for checking consistency of only OPRA constraints  (this file used for OPRA experiments)


disjunctive.lp :  ASP subprogram for disjunctive HOPA  constraints

infer.lp :  ASP subprogram for inference

inconsist.lp :  ASP subprogram for explaining inconsistency

presumed.lp :  ASP subprogram for presumed constraints

tangent.lp :  includes tangent values of angles in the range [0,90] degree

sqrtX.lp : includes square root of integers in the range [0,X]. We need this because ASP solver clingo itself does not compute square root. 

gridS.lp :  specifies a grid of size SxS,  and corresponding distance intervals (upper and lower bound) for qualitative distance relations

degZ.lp :  specifies angle resolution as Z degrees  and corresponding angle values (for choosing orientation of objects)



SCENARIOS:

The folder "Scenario" includes the ASP encoding of the two scenarios in the paper.


To obtain the output for Marine Navigation scenario, run

./clingo  ./ASP/hopa_main.lp  ./ASP/presumed.lp  ./ASP/infer.lp  ./ASP/tangent.lp  ./ASP/sqrt50.lp  ./Scenario/marine.lp


To obtain the output for Robotic Perception scenario, run

./clingo  ./ASP/hopa_main.lp  ./ASP/inconsist.lp  ./ASP/tangent.lp  ./ASP/sqrt50.lp  ./Scenario/robot.lp



HOPA BENCHMARK TEST INSTANCES:

The folder "HOPA_Instances"  involves basic and disjunctive (in)consistent HOPA benchmark test instances used in the experiments (100 samples for each). 

"nX_cY_gridZ_basic_sat.lp" is the consistent problem instance with basic HOPA constraints, X objects, Y constraints and grid size Z

"nX_cY_gridZ_basic_unsat.lp" is the inconsistent problem instance with basic HOPA constraints, X objects, Y constraints and grid size Z


"nX_cY_gridZ_disjDxF_sat.lp" is the consistent disjunctive HOPA problem instance with D disjunctive constraints each having F disjuncts

"nX_cY_gridZ_disjDxF_unsat.lp" is the inconsistent disjunctive HOPA problem instance with D disjunctive constraints each having F disjuncts



EXAMPLES: To test a basic consistent HOPA problem instance, open a terminal in the main folder and run

./clingo  ./ASP/hopa_main.lp  ./ASP/grid50.lp  ./ASP/deg3.lp   ./ASP/sqrt50.lp  ./ASP/tangent.lp   ./HOPA_Instances/Consistent/Sample1/Basic/n40_c40_grid50_basic_sat.lp

(Can make test with other angular resolution e.g. deg2.lp for resolution 2 degrees)


To test a basic inconsistent HOPA problem instance, open a terminal in the main folder and run

./clingo  ./ASP/hopa_main.lp  ./ASP/grid50.lp  ./ASP/deg3.lp   ./ASP/sqrt50.lp  ./ASP/tangent.lp   ./HOPA_Instances/Inconsistent/Sample1/Basic/n40_c40_grid50_basic_unsat.lp

(Can make test with other angular resolution e.g. deg2.lp for resolution 2 degrees)


To test a disjunctive consistent HOPA problem instance, open a terminal in the main folder and run

./clingo  ./ASP/hopa_main.lp  ./ASP/grid50.lp  ./ASP/deg3.lp   ./ASP/sqrt50.lp  ./ASP/tangent.lp   ./HOPA_Instances/Consistent/Sample1/Disjunctive/n40_c40_grid50_disj8x4_sat.lp

(Can make test with other angular resolution e.g. deg2.lp for resolution 2 degrees)



To test a disjunctive inconsistent HOPA problem instance, open a terminal in the main folder and run

./clingo  ./ASP/hopa_main.lp  ./ASP/grid50.lp  ./ASP/deg3.lp   ./ASP/sqrt50.lp  ./ASP/tangent.lp   ./HOPA_Instances/Inconsistent/Sample1/Disjunctive/n40_c40_grid50_disj8x4_unsat.lp

(Can make test with other angular resolution e.g. deg2.lp for resolution 2 degrees)



OPRA BENCHMARK TEST INSTANCES:

The folder "OPRA_Instances"  involves basic (in)consistent OPRA benchmark test instances used in the experiments (100 samples for each). 

"nX_cY_gZ_basic_sat.lp" is the consistent problem instance with basic OPRA constraints, X objects, Y constraints and granularity Z

"nX_cY_gZ_basic_unsat.lp" is the inconsistent problem instance with basic OPRA constraints, X objects, Y constraints and granularity Z

Note: "opra_consistency.lp" is the simplified version of the main ASP program HOPA, for checking consistency of only OPRA constraints. 



EXAMPLES: To test a basic consistent OPRA problem instance, open a terminal in the main folder and run

./clingo  ./ASP/opra_consistency.lp  ./ASP/grid50.lp  ./ASP/deg3.lp  ./ASP/tangent.lp   ./OPRA_Instances/Consistent/Sample1/Basic/n20_c20_g2_basic_sat.lp

(Can make test with other grid size or angular resolution; e.g. grid100.lp for grid size 100, and deg2.lp for resolution 2 degrees)


To test a basic inconsistent OPRA problem instance, open a terminal in the main folder and run

./clingo  ./ASP/opra_consistency.lp  ./ASP/grid50.lp  ./ASP/deg3.lp  ./ASP/tangent.lp   ./OPRA_Instances/Inconsistent/Sample1/Basic/n20_c20_g2_basic_unsat.lp

(Can make test with other grid size or angular resolution; e.g. grid100.lp for grid size 100, and deg2.lp for resolution 2 degrees)


