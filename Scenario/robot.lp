
% Robotic Example


% NOTATION:  
% diff_opra(N,K,L,I,J,G):   object K is on segment I wrt L, object L is on segment J wrt K
% same_opra(N,K,L,I,G):     object L is on segment I wrt K
% one_sided_diff_opra(N,K,L,I,G):  object K is on segment I wrt L



% Grid size

xcoord(1..15).
ycoord(1..10).



% Angle and resolution

angle_res(3).

angle_range(0..119).



% Points:  1: r1,  2: r2,  3: b1,  4:b2,  5: c 

object(1..11).  

point(1,1).  point(2,1).  point(3,1).  point(4,1).  point(5,2).  point(6,2).  point(7,2).  point(8,2).  point(9,3).  point(10,4).  point(11,5).



% Constraint ID

constraint_id(1..19).

mandatory(1..17).




% Object's orientations are stated as vector

direction(1,1,0).

direction(2,1,1,2).  direction(3,2,1,3).  direction(4,3,1,4).  direction(5,4,1,5).  direction(6,5,1,2).  direction(7,6,2,3).  direction(8,7,2,4).  direction(9,8,2,5).  direction(10,11,5,1).




% Initial perception at point R1

same_opra(11,1,4,1,4).

same_opra(12,1,2,1,4).

same_opra(13,1,3,15,4).



% Perception at point R2

diff_opra(14,1,5,8,0,4).

same_opra(15,5,8,4,4).

same_opra(16,5,6,1,4).

same_opra(17,5,7,13,4).



% Check which object is behind:

% disjrelation(N,D,3,K,L,I,G).

disj_index(18,1..3).
disjrelation(18,1,2,9,11,7,4).
disjrelation(18,2,2,9,11,8,4).
disjrelation(18,3,2,9,11,9,4).

disj_index(19,1..3).
disjrelation(19,1,2,10,11,7,4).
disjrelation(19,2,2,10,11,8,4).
disjrelation(19,3,2,10,11,9,4).


