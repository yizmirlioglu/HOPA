

% NOTATION:  
% diff_opra(N,K,L,I,J,G):   object K is on segment I wrt L, object L is on segment J wrt K
% same_opra(N,K,L,I,G):     object L is on segment I wrt K
% one_sided_diff_opra(N,K,L,I,G):  object K is on segment I wrt L

% direction(K,A,B) means orientation of object K is from point A to B
% direction(K,T) means orientation of object K is known and equal to T
% both angle direction and angle vector can be specified together

% Other Quantitative HOPA constraints are represented by psi_precise, psi_range, delta_precise, delta_range, location, distance, distance_range atoms.

% Qualitative distance constraint is represented by dist_rel atom.


numobj(N) :- N=#count{K:object(K)}.



% Disjunctive constraints
% Choose one basic relation from every disjunctive constraint

% 1: same opra, 2: one-sided diff opra 3: diff opra 4: qualit dist, 5:orient exact, 6:orient range, 7:orient vector, 8.psi exact, 9:psi range, 10:delta exact, 11. delt range, 12:exact location  13: distance range, 14: exact distance

existdisj(N) :- not drop(N), disj_index(N,D). 

1{chosen(N,D): disj_index(N,D)}1 :- not drop(N), existdisj(N).


same_opra(N,K,L,I,G) :- chosen(N,D), disjrelation(N,D,1,K,L,I,G).

one_sided_diff_opra(N,K,L,I,G) :- chosen(N,D), disjrelation(N,D,2,K,L,I,G).

diff_opra(N,K,L,I,J,G) :- chosen(N,D), disjrelation(N,D,3,K,L,I,J,G).

dist_rel(N,K,L,U) :- chosen(N,D), disjrelation(N,D,4,K,L,U).

direction(N,K,A) :- chosen(N,D), disjrelation(N,D,5,K,A).

direction_range(N,K,A1,A2) :- chosen(N,D), disjrelation(N,D,6,K,A1,A2).

direction(N,K,P,R) :- chosen(N,D), disjrelation(N,D,7,K,P,R).

psi_precise(N,K,L,A) :- chosen(N,D), disjrelation(N,D,8,K,L,A).

psi_range(N,K,L,A1,A2) :- chosen(N,D), disjrelation(N,D,9,K,L,A1,A2).

delta_precise(N,K,L,A) :- chosen(N,D), disjrelation(N,D,10,K,L,A).

delta_range(N,K,L,A1,A2) :- chosen(N,D), disjrelation(N,D,11,K,L,A1,A2).

location(N,K,X,Y) :- chosen(N,D), disjrelation(N,D,12,K,X,Y).

distance_range(N,K,L,U1,U2) :- chosen(N,D), disjrelation(N,D,13,K,L,U1,U2).

distance(N,K,L,U) :- chosen(N,D), disjrelation(N,D,14,K,L,U).






% Transform differential OPRA constraints to one-sided differential OPRA constraints

single_diff_opra(N,K,L,I,G) :- not drop(N), diff_opra(N,K,L,I,J,G).

single_diff_opra(N,L,K,J,G) :- not drop(N), diff_opra(N,K,L,I,J,G).


single_diff_opra(N,K,L,I,G) :- not drop(N), one_sided_diff_opra(N,K,L,I,G).






% Check whether same and differential OPRA constraints are specified for same pair (for detecting inconsistency faster)

:- not drop(N), single_diff_opra(M,K,L,I,G1), same_opra(N,L,K,J,G2).

:- not drop(N), single_diff_opra(M,K,L,I,G1), same_opra(N,K,L,J,G2).

:- I!=J, single_diff_opra(M,K,L,I,G), single_diff_opra(N,K,L,J,G).

:- I!=0, J!=0, I+J!=4*G, not drop(M), not drop(N), same_opra(M,K,L,I,G), same_opra(N,L,K,J,G).




% Instantiate objects which appear in some OPRA constraint

instantiate(P) :- not drop(N), point(K,P), diff_opra(N,K,L,I,J,G).
instantiate(R) :- not drop(N), point(L,R), diff_opra(N,K,L,I,J,G).

instantiate(P) :- not drop(N), point(K,P), same_opra(N,K,L,I,G).
instantiate(R) :- not drop(N), point(L,R), same_opra(N,K,L,I,G).

instantiate(P) :- not drop(N), point(K,P), one_sided_diff_opra(N,K,L,I,G).
instantiate(R) :- not drop(N), point(L,R), one_sided_diff_opra(N,K,L,I,G).

instantiate(P) :- not drop(N), point(K,P), location(N,K,X,Y).

instantiate(P) :- not drop(N), point(K,P), direction(N,K,A).
instantiate(P) :- not drop(N), direction(N,K,P,R).
instantiate(R) :- not drop(N), direction(N,K,P,R).
instantiate(P) :- not drop(N), point(K,P), direction_range(N,K,A1,A2).

instantiate(P) :- not drop(N), point(K,P), psi_precise(N,K,L,A).
instantiate(R) :- not drop(N), point(L,R), psi_precise(N,K,L,A).
instantiate(P) :- not drop(N), point(K,P), psi_range(N,K,L,A1,A2).
instantiate(R) :- not drop(N), point(L,R), psi_range(N,K,L,A1,A2).

instantiate(P) :- not drop(N), point(K,P), delta_precise(N,K,L,A).
instantiate(R) :- not drop(N), point(L,R), delta_precise(N,K,L,A).
instantiate(P) :- not drop(N), point(K,P), delta_range(N,K,L,A1,A2).
instantiate(R) :- not drop(N), point(L,R), delta_range(N,K,L,A1,A2).

instantiate(P) :- not drop(N), point(K,P), dist_rel(N,K,L,U).
instantiate(R) :- not drop(N), point(L,R), dist_rel(N,K,L,U).

instantiate(P) :- not drop(N), point(K,P), distance(N,K,L,D).
instantiate(R) :- not drop(N), point(L,R), distance(N,K,L,D).
instantiate(P) :- not drop(N), point(K,P), distance_range(N,K,L,D1,D2).
instantiate(R) :- not drop(N), point(L,R), distance_range(N,K,L,D1,D2).




% Objects to generate their orientation

generate_orient(K) :- not drop(N), diff_opra(N,K,L,I,J,G).

generate_orient(L) :- not drop(N), diff_opra(N,K,L,I,J,G).

generate_orient(L) :- not drop(N), one_sided_diff_opra(N,K,L,I,G).

generate_orient(K) :- not drop(N), same_opra(N,K,L,I,G).

generate_orient(L) :- not drop(N), same_opra(N,K,L,I,G).


generate_orient(L) :- not drop(N), psi_precise(N,K,L,A).

generate_orient(L) :- not drop(N), psi_range(N,K,L,A1,A2).

generate_orient(K) :- not drop(N), delta_precise(N,K,L,A).
generate_orient(L) :- not drop(N), delta_precise(N,K,L,A).

generate_orient(K) :- not drop(N), delta_range(N,K,L,A1,A2).
generate_orient(L) :- not drop(N), delta_range(N,K,L,A1,A2).


generate_orient(K) :- not drop(N), direction(N,K,A).

generate_orient(K) :- not drop(N), direction(N,K,P,R).

generate_orient(K) :- not drop(N), direction_range(N,K,A1,A2).






% Objects whose Locations are already known or can be deduced

loc_known(P) :- not drop(N), point(K,P), location(N,K,X,Y).

loc_known(P) :- P<R, not drop(N), point(K,P), point(L,R), same_opra(N,K,L,I,G).
loc_known(R) :- P>R, not drop(N), point(K,P), point(L,R), same_opra(N,K,L,I,G).

loc_known(P) :- P<R, not drop(N), point(K,P), point(L,R), delta_precise(N,K,L,A).
loc_known(R) :- P>R, not drop(N), point(K,P), point(L,R), delta_precise(N,K,L,A).

loc_known(P) :- P<R, not drop(N), point(K,P), point(L,R), delta_range(N,K,L,A1,A2).
loc_known(R) :- P>R, not drop(N), point(K,P), point(L,R), delta_range(N,K,L,A1,A2).





% Objects whose orientation information is specified by user or can be implicitly computed

orient_specified(K) :- T>=0, not drop(N), direction(N,K,T), object(K).

orient_specified(K) :- not drop(N), direction(N,K,P,R), object(K).

orient_specified(L) :- I\2=0, single_diff_opra(N,K,L,I,G).

orient_specified(K) :- I\2=0, K<L, not drop(N), same_opra(N,K,L,I,G).
orient_specified(L) :- I\2=0, K>L, not drop(N), same_opra(N,K,L,I,G).

orient_specified(K) :- K<L, not drop(N), delta_precise(N,K,L,A).
orient_specified(L) :- K>L, not drop(N), delta_precise(N,K,L,A).

orient_specified(L) :- not drop(N), psi_precise(N,K,L,T).



% Exogenous numerical distance information 

dist_known(P,R) :- not drop(N), point(K,P), point(L,R), distance(N,K,L,D).





% Choose locations

1{xloc(P,X): xcoord(X)}1 :- not loc_known(P), instantiate(P).  

1{yloc(P,Y): ycoord(Y)}1 :- not loc_known(P), instantiate(P). 




% Forbidden locations and obstacles

:- invalid_loc(X,Y), xloc(P,X), yloc(P,Y), instantiate(P). 





% Locations of objects that are specified by the user or implicitly deduced

xloc(P,X) :- not drop(N), point(K,P), location(N,K,X,Y).
yloc(P,Y) :- not drop(N), point(K,P), location(N,K,X,Y).

xloc(P,X) :- P<R, xloc(R,X), not drop(N), point(K,P), point(L,R), same_opra(N,K,L,I,G).
xloc(R,X) :- P>R, xloc(P,X), not drop(N), point(K,P), point(L,R), same_opra(N,K,L,I,G).

yloc(P,Y) :- P<R, yloc(R,Y), not drop(N), point(K,P), point(L,R), same_opra(N,K,L,I,G).
yloc(R,Y) :- P>R, yloc(P,Y), not drop(N), point(K,P), point(L,R), same_opra(N,K,L,I,G).


xloc(P,X) :- P<R, xloc(R,X), not drop(N), point(K,P), point(L,R), delta_precise(N,K,L,A).
xloc(R,X) :- P>R, xloc(P,X), not drop(N), point(K,P), point(L,R), delta_precise(N,K,L,A).

yloc(P,Y) :- P<R, yloc(R,Y), not drop(N), point(K,P), point(L,R), delta_precise(N,K,L,A).
yloc(R,Y) :- P>R, yloc(P,Y), not drop(N), point(K,P), point(L,R), delta_precise(N,K,L,A).


xloc(P,X) :- P<R, xloc(R,X), not drop(N), point(K,P), point(L,R), delta_range(N,K,L,A1,A2).
xloc(R,X) :- P>R, xloc(P,X), not drop(N), point(K,P), point(L,R), delta_range(N,K,L,A1,A2).

yloc(P,Y) :- P<R, yloc(R,Y), not drop(N), point(K,P), point(L,R), delta_range(N,K,L,A1,A2).
yloc(R,Y) :- P>R, yloc(P,Y), not drop(N), point(K,P), point(L,R), delta_range(N,K,L,A1,A2).





% Compute vector between objects

compute_vector(P,R) :- P>R, not drop(N), point(K,P), point(L,R), diff_opra(N,K,L,I,J,G).
compute_vector(R,P) :- P<R, not drop(N), point(K,P), point(L,R), diff_opra(N,K,L,I,J,G).

compute_reverse(P,R) :- P>R, not drop(N), point(K,P), point(L,R), diff_opra(N,K,L,I,J,G).
compute_reverse(R,P) :- P<R, not drop(N), point(K,P), point(L,R), diff_opra(N,K,L,I,J,G).

compute_vector(P,R) :- not drop(N), point(K,P), point(L,R), one_sided_diff_opra(N,K,L,I,G).


compute_vector(P,R) :- not drop(N), point(K,P), point(L,R), psi_precise(N,K,L,A).

compute_vector(P,R) :- not drop(N), point(K,P), point(L,R), psi_range(N,K,L,A1,A2).





% Compute vector corresponding to direction

compute_vector(R,P) :- not drop(N), direction(N,L,P,R).

compute_vector(R,P) :- not drop(N), direction(N,L,P,R), single_diff_opra(M,K,L,I,G).

compute_vector(R,P) :- not drop(N), not drop(M), direction(N,K,P,R), same_opra(M,K,L,I,G).

compute_vector(R,P) :- not drop(N), not drop(M), direction(N,L,P,R), same_opra(M,K,L,I,G).





% Compute vector between objects

xdist(P,R,X1-X2) :- xloc(P,X1), xloc(R,X2), compute_vector(P,R).
ydist(P,R,Y1-Y2) :- yloc(P,Y1), yloc(R,Y2), compute_vector(P,R).

xdist_abs(P,R,|DX|) :- xdist(P,R,DX), compute_vector(P,R).
ydist_abs(P,R,|DY|) :- ydist(P,R,DY), compute_vector(P,R).

xdist_positive(P,R) :- DX>=0, xdist(P,R,DX), compute_vector(P,R).
ydist_positive(P,R) :- DY>=0, ydist(P,R,DY), compute_vector(P,R).




% Compute Angle of the Vector

qangle(P,R,0) :- ydist(P,R,0), compute_vector(P,R).

qangle(P,R,90) :- xdist(P,R,0), compute_vector(P,R).


tan_abs(P,R,C*DY/DX) :- DX>0, DY>0, tan_coeff(C), xdist_abs(P,R,DX), ydist_abs(P,R,DY), compute_vector(P,R).

qangle(P,R,A) :- T>=R1, T<=R2, tan_abs(P,R,T), tan_range(A,R1,R2), compute_vector(P,R).



varphi(P,R,A) :- qangle(P,R,A), xdist_positive(P,R), ydist_positive(P,R), compute_vector(P,R).

varphi(P,R,180-A) :- qangle(P,R,A), not xdist_positive(P,R), ydist_positive(P,R), compute_vector(P,R).

varphi(P,R,360-A) :- qangle(P,R,A), xdist_positive(P,R), not ydist_positive(P,R), compute_vector(P,R).

varphi(P,R,180+A) :- qangle(P,R,A), not xdist_positive(P,R), not ydist_positive(P,R), compute_vector(P,R).


varphi(R,P,(A+180)\360) :- varphi(P,R,A), compute_reverse(P,R).






% Choose orientation and quadrant for point objects

1{degree(K,D): angle_range(D)}1 :- not orient_specified(K), generate_orient(K), object(K).     % direction(K,-1)

orient(K,D*R) :- degree(K,D), angle_res(R), not orient_specified(K), generate_orient(K), object(K).      % direction(K,-1)

orient(K,T) :- T>=0, not drop(N), direction(N,K,T), object(K).

orient(K,A) :- not drop(N), varphi(R,P,A), direction(N,K,P,R), object(K).

orient(L,(360+A-I*90/G)\360) :- I\2=0, varphi(P,R,A), point(K,P), point(L,R), single_diff_opra(N,K,L,I,G).

orient(L,(A-T+360)\360) :- not drop(N), varphi(P,R,A), point(K,P), point(L,R), psi_precise(N,K,L,T).

orient(K,(360+T-I*90/G)\360) :- I\2=0, K<L, orient(L,T), not drop(N), same_opra(N,K,L,I,G).
orient(L,(360+T+I*90/G)\360) :- I\2=0, K>L, orient(K,T), not drop(N), same_opra(N,K,L,I,G).

orient(K,(T-A+360)\360) :- K<L, orient(L,T), not drop(N), delta_precise(N,K,L,A).
orient(L,(T+A+360)\360) :- K>L, orient(K,T), not drop(N), delta_precise(N,K,L,A).


:- I\2=0, T>=0, (A-T+360)\360!=I*90/G, not drop(N), direction(N,L,T), varphi(P,R,A), point(K,P), point(L,R), single_diff_opra(M,K,L,I,G).

:- I\2=0, (A-T+360)\360!=I*90/G, varphi(Z,Y,T), not drop(N), direction(N,L,Y,Z), varphi(P,R,A), point(K,P), point(L,R), single_diff_opra(M,K,L,I,G).

:- T>=0, T!=A, varphi(R,P,A), not drop(N), not drop(M), direction(N,K,T), direction(M,K,P,R), object(K).

:- T1!=T2, orient(K,T1), orient(K,T2).





% Differential OPRA constraints

% First impose that locations are different

:- xdist(P,R,0), ydist(P,R,0), point(K,P), point(L,R), single_diff_opra(N,K,L,I,G).




% Angular differential OPRA constraints

% impose orientation when bounds normal

:- T<=(360+A-(I+1)*90/G)\360, I\2=1, (360+A-(I-1)*90/G)\360>(360+A-(I+1)*90/G)\360, orient(L,T), varphi(P,R,A), point(K,P), point(L,R), single_diff_opra(N,K,L,I,G).

:- T>=(360+A-(I-1)*90/G)\360, I\2=1, (360+A-(I-1)*90/G)\360>(360+A-(I+1)*90/G)\360, orient(L,T), varphi(P,R,A), point(K,P), point(L,R), single_diff_opra(N,K,L,I,G).


% impose orientation when bounds reverse

:- T<=(360+A-(I+1)*90/G)\360, T>=(360+A-(I-1)*90/G)\360, I\2=1, (360+A-(I-1)*90/G)\360<(360+A-(I+1)*90/G)\360, orient(L,T), varphi(P,R,A), point(K,P), point(L,R), single_diff_opra(N,K,L,I,G).




% Check orientation of objects related by Same OPRA constraints

:- I\2=1, (T-B+360)\360<=(I-1)*90/G, orient(K,B), orient(L,T), not drop(N), same_opra(N,K,L,I,G).

:- I\2=1, (T-B+360)\360>=(I+1)*90/G, orient(K,B), orient(L,T), not drop(N), same_opra(N,K,L,I,G).




% Quantitative Constraints

:- T<A1, not drop(N), orient(K,T), direction_range(N,K,A1,A2), object(K).

:- T>A2, not drop(N), orient(K,T), direction_range(N,K,A1,A2), object(K).

:- xdist(P,R,0), ydist(P,R,0), point(K,P), point(L,R), not drop(N), psi_precise(N,K,L,T).

:- xdist(P,R,0), ydist(P,R,0), point(K,P), point(L,R), not drop(N), psi_range(N,K,L,A1,A2).

:- T>(A-A1+360)\360, (A-A1+360)\360>(A-A2+360)\360, orient(L,T), varphi(P,R,A), point(K,P), point(L,R), not drop(N), psi_range(N,K,L,A1,A2).

:- T<(A-A2+360)\360, (A-A1+360)\360>(A-A2+360)\360, orient(L,T), varphi(P,R,A), point(K,P), point(L,R), not drop(N), psi_range(N,K,L,A1,A2).

:- T>(A-A1+360)\360, T<(A-A2+360)\360, (A-A1+360)\360<(A-A2+360)\360, orient(L,T), varphi(P,R,A), point(K,P), point(L,R), not drop(N), psi_range(N,K,L,A1,A2).

:- (T-B+360)\360<A1, orient(K,B), orient(L,T), not drop(N), delta_range(N,K,L,A1,A2).

:- (T-B+360)\360>A2, orient(K,B), orient(L,T), not drop(N), delta_range(N,K,L,A1,A2).






% DISTANCE  CONSTRAINTS

compute_loc_diff(P,R) :- not compute_vector(P,R), point(K,P), point(L,R), not drop(N), dist_rel(N,K,L,U).  % not dist_given(P,R), 

compute_loc_diff(P,R) :- not compute_vector(P,R), point(K,P), point(L,R), not drop(N), distance_range(N,K,L,D1,D2).


compute_dist(P,R) :- not dist_known(P,R), point(K,P), point(L,R), not drop(N), dist_rel(N,K,L,U). 

compute_dist(P,R) :- not dist_known(P,R), point(K,P), point(L,R), not drop(N), distance_range(N,K,L,D1,D2). 


compute_loc_diff(P,R) :- not compute_vector(P,R), point(K,P), point(L,R), not drop(N), distance(N,K,L,D).





% Compute Numerical Distance

xdist(P,R,X1-X2) :- xloc(P,X1), xloc(R,X2), compute_loc_diff(P,R).
ydist(P,R,Y1-Y2) :- yloc(P,Y1), yloc(R,Y2), compute_loc_diff(P,R).

ndist(P,R,DX*DX+DY*DY) :- xdist(P,R,DX), ydist(P,R,DY), compute_dist(P,R).

ndist(P,R,D) :- point(K,P), point(L,R), not drop(N), distance(N,K,L,D).




% Check Qualitative Distance Constraints

% Symmetry

:- U1!=U2, not drop(M), not drop(N), dist_rel(M,K,L,U1), dist_rel(N,K,L,U2).

:- U1!=U2, not drop(M), not drop(N), dist_rel(M,K,L,U1), dist_rel(N,L,K,U2).



% Check interval for qualitative distance constraints

:- D<D1, ndist(P,R,D), lower_bound(U,D1), point(K,P), point(L,R), not drop(N), dist_rel(N,K,L,U).

:- D>D1, ndist(P,R,D), upper_bound(U,D1), point(K,P), point(L,R), not drop(N), dist_rel(N,K,L,U).




% Check Quantitative Distance Constraints

ydist_known(P,R,Y) :- sqroot(D-X*X,Y), xdist(P,R,X), point(K,P), point(L,R), not drop(N), distance(N,K,L,D).

ydist_known(P,R) :- ydist_known(P,R,Y).

:- not ydist_known(P,R), point(K,P), point(L,R), not drop(N), distance(N,K,L,D).

:- Y1!=Y2-DY, Y1!=Y2+DY, ydist_known(P,R,DY), yloc(P,Y1), yloc(R,Y2), point(K,P), point(L,R), not drop(N), distance(N,K,L,D).

:- D<D1, ndist(P,R,D), point(K,P), point(L,R), not drop(N), distance_range(N,K,L,D1,D2).

:- D>D2, ndist(P,R,D), point(K,P), point(L,R), not drop(N), distance_range(N,K,L,D1,D2).



#show xloc/2.
#show yloc/2.
#show orient/2.
#show chosen/2.
