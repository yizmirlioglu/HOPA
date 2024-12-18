
% Marine Navigation Scenario


% NOTATION:  
% diff_opra(N,K,L,I,J,G):   object K is on segment I wrt L, object L is on segment J wrt K
% same_opra(N,K,L,I,G):     object L is on segment I wrt K
% one_sided_diff_opra(N,K,L,I,G):  object K is on segment I wrt L

% direction(K,T), T>=0 means orientation of opoint K is known and T


% Size of the Grid

xcoord(1..40).
ycoord(1..30).


% Angle and resolution

angle_res(3).

angle_range(0..119).



% 1: Port station, 2..7: ships,  8: Kumburun point

object(1..8).  

point(1,1).  point(2,2).  point(3,3).  point(4,4).  point(5,5).  point(6,6).  point(7,7).   point(8,8). 


% Infer

toinfer_orient(3).

toinfer_opra(3,4,4).



% Constraint ID

constraint_id(1..13).

presumed_constraint_id(14).



% Forbidden locations (land area)

invalid_loc(37..40,1..30).  
invalid_loc(33..36,1..6).  
invalid_loc(34..36,9..18). 
invalid_loc(34..36,23..27). 
invalid_loc(28..31,21..24). 
invalid_loc(30..31,18..21). 



% Constraints

location(1,1,33,16).

location(2,2,32,6).

location(3,8,33,24).

psi_precise(4,2,4,213).

dist_rel(5,2,4,near). 

one_sided_diff_opra(6,4,5,7,2).

location(7,5,18,11).

distance_range(8,6,1,400,625).

disj_index(9,1..5).
disjrelation(9,1,2,6,3,6,4).
disjrelation(9,2,2,6,3,7,4).
disjrelation(9,3,2,6,3,8,4).
disjrelation(9,4,2,6,3,9,4).
disjrelation(9,5,2,6,3,10,4).

disj_index(10,1..3).
disjrelation(10,1,4,6,3,medium).
disjrelation(10,2,4,6,3,far).
disjrelation(10,3,4,6,3,veryfar).

one_sided_diff_opra(11,6,7,6,6).

dist_rel(12,6,7,verynear).

psi_range(13,5,3,70,100).




% presumed constraint

presumed_direction(14,6,6,8).




% Distance Relations Bounds:  Interval for each distance relation  (square of distance)

lower_bound(on,0).  lower_bound(verynear,1).  lower_bound(near,26).  lower_bound(medium,101).  lower_bound(far,401).  lower_bound(veryfar,1226).

upper_bound(on,0).  upper_bound(verynear,25).  upper_bound(near,100).  upper_bound(medium,400).  upper_bound(far,1225).


