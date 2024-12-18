

% This is the simplified version of the main HOPA ASP program which is used for checking consistency of only OPRA constraints. 
% Since object locations are unknown, point identifier is unnecessary and removed.



% NOTATION:  
% diff_opra(N,K,L,I,J,G):   object K is on segment I wrt L, object L is on segment J wrt K
% same_opra(N,K,L,I,G):     object L is on segment I wrt K
% one_sided_diff_opra(N,K,L,I,G):  object K is on segment I wrt L


numobj(N) :- N=#count{K:object(K)}.




% Transform differential OPRA constraints to one-sided differential OPRA constraints

single_diff_opra(N,K,L,I,G) :- not drop(N), diff_opra(N,K,L,I,J,G).

single_diff_opra(N,L,K,J,G) :- not drop(N), diff_opra(N,K,L,I,J,G).


single_diff_opra(N,K,L,I,G) :- not drop(N), one_sided_diff_opra(N,K,L,I,G).





% Impose symmetry for consistency

:- not drop(N), single_diff_opra(M,K,L,I,G1), same_opra(N,L,K,J,G2).

:- not drop(N), single_diff_opra(M,K,L,I,G1), same_opra(N,K,L,J,G2).

:- I!=J, single_diff_opra(M,K,L,I,G), single_diff_opra(N,K,L,J,G).

:- I!=0, J!=0, I+J!=4*G, not drop(M), not drop(N), same_opra(M,K,L,I,G), same_opra(N,L,K,J,G).





% Instantiate objects which appear in some OPRA constraint

instantiate(K) :- not drop(N), diff_opra(N,K,L,I,J,G).
instantiate(L) :- not drop(N), diff_opra(N,K,L,I,J,G).

instantiate(K) :- not drop(N), same_opra(N,K,L,I,G).
instantiate(L) :- not drop(N), same_opra(N,K,L,I,G).

instantiate(K) :- not drop(N), one_sided_diff_opra(N,K,L,I,G).
instantiate(L) :- not drop(N), one_sided_diff_opra(N,K,L,I,G).





% Compute the vector between objects

compute_vector(K,L) :- K>L, not drop(N), diff_opra(N,K,L,I,J,G).
compute_vector(L,K) :- K<L, not drop(N), diff_opra(N,K,L,I,J,G).

compute_reverse(K,L) :- K>L, not drop(N), diff_opra(N,K,L,I,J,G).
compute_reverse(L,K) :- K<L, not drop(N), diff_opra(N,K,L,I,J,G).

compute_vector(K,L) :- not drop(N), one_sided_diff_opra(N,K,L,I,G).





% Objects whose Locations are already known or can be deduced

loc_known(K) :- K<L, not drop(N), same_opra(N,K,L,I,G).
loc_known(L) :- K>L, not drop(N), same_opra(N,K,L,I,G).




% Objects whose orientation information is specified by user or can be implicitly computed

orient_specified(L) :- I\2=0, single_diff_opra(N,K,L,I,G).

orient_specified(K) :- I\2=0, K<L, not drop(N), same_opra(N,K,L,I,G).
orient_specified(L) :- I\2=0, K>L, not drop(N), same_opra(N,K,L,I,G).





% Choose locations

1{xloc(K,X): xcoord(X)}1 :- not loc_known(K), instantiate(K).

1{yloc(K,Y): ycoord(Y)}1 :- not loc_known(K), instantiate(K).





% Locations for Same OPRA relations

xloc(K,X) :- K<L, xloc(L,X), not drop(N), same_opra(N,K,L,I,G).
xloc(L,X) :- K>L, xloc(K,X), not drop(N), same_opra(N,K,L,I,G).

yloc(K,Y) :- K<L, yloc(L,Y), not drop(N), same_opra(N,K,L,I,G).
yloc(L,Y) :- K>L, yloc(K,Y), not drop(N), same_opra(N,K,L,I,G).





% Choose orientation and quadrant for objects

1{degree(K,D): angle_range(D)}1 :- not orient_specified(K), instantiate(K).

orient(K,D*R) :- degree(K,D), angle_res(R), not orient_specified(K), instantiate(K).



% Find alreadt known or deduced orientation 

orient(L,(360+A-I*90/G)\360) :- I\2=0, varphi(K,L,A), single_diff_opra(N,K,L,I,G).

orient(K,(360+T-I*90/G)\360) :- I\2=0, K<L, orient(L,T), not drop(N), same_opra(N,K,L,I,G).
orient(L,(360+T+I*90/G)\360) :- I\2=0, K>L, orient(K,T), not drop(N), same_opra(N,K,L,I,G).

:- A1!=A2, orient(L,A1), orient(L,A2).




% Compute vector between objects

xdist(K,L,X1-X2) :- xloc(K,X1), xloc(L,X2), compute_vector(K,L).
ydist(K,L,Y1-Y2) :- yloc(K,Y1), yloc(L,Y2), compute_vector(K,L).


xdist_abs(K,L,|DX|) :- xdist(K,L,DX), compute_vector(K,L).
ydist_abs(K,L,|DY|) :- ydist(K,L,DY), compute_vector(K,L).

xdist_positive(K,L) :- DX>=0, xdist(K,L,DX), compute_vector(K,L).
ydist_positive(K,L) :- DY>=0, ydist(K,L,DY), compute_vector(K,L).




% Compute Angle of the Vector


qangle(K,L,0) :- ydist(K,L,0), compute_vector(K,L).

qangle(K,L,90) :- xdist(K,L,0), compute_vector(K,L).

tan_abs(K,L,C*DY/DX) :- DX>0, DY>0, tan_coeff(C), xdist_abs(K,L,DX), ydist_abs(K,L,DY), compute_vector(K,L).

qangle(K,L,A) :- T>=R1, T<=R2, tan_abs(K,L,T), tan_range(A,R1,R2), compute_vector(K,L).



varphi(K,L,A) :- qangle(K,L,A), xdist_positive(K,L), ydist_positive(K,L), compute_vector(K,L).

varphi(K,L,360-A) :- qangle(K,L,A), xdist_positive(K,L), not ydist_positive(K,L), compute_vector(K,L).

varphi(K,L,180-A) :- qangle(K,L,A), not xdist_positive(K,L), ydist_positive(K,L), compute_vector(K,L).

varphi(K,L,180+A) :- qangle(K,L,A), not xdist_positive(K,L), not ydist_positive(K,L), compute_vector(K,L).


varphi(L,K,(A+180)\360) :- varphi(K,L,A), compute_reverse(K,L).





% Differential OPRA constraints

% First, impose locations cannot be the same

:- xdist(K,L,0), ydist(K,L,0), single_diff_opra(N,K,L,I,G).

%:- K>L, xdist(K,L,0), ydist(K,L,0), single_diff_opra(N,K,L,I,G).
%:- K<L, xdist(L,K,0), ydist(L,K,0), single_diff_opra(N,K,L,I,G).




% Find range of orientation for each object for Non-linear diff OPRA constraints

% impose orientation when bounds normal

:- P<=(360+A-(I+1)*90/G)\360, I\2=1, (360+A-(I-1)*90/G)\360>(360+A-(I+1)*90/G)\360, orient(L,P), varphi(K,L,A), single_diff_opra(N,K,L,I,G).

:- P>=(360+A-(I-1)*90/G)\360, I\2=1, (360+A-(I-1)*90/G)\360>(360+A-(I+1)*90/G)\360, orient(L,P), varphi(K,L,A), single_diff_opra(N,K,L,I,G).



% impose orientation when bounds reverse

:- P<=(360+A-(I+1)*90/G)\360, P>=(360+A-(I-1)*90/G)\360, I\2=1, (360+A-(I-1)*90/G)\360<(360+A-(I+1)*90/G)\360, orient(L,P), varphi(K,L,A), single_diff_opra(N,K,L,I,G).




% Check orientation of objects related by same OPRA constraints

:- I\2=1, (T-P+360)\360<=(I-1)*90/G, orient(K,P), orient(L,T), not drop(N), same_opra(N,K,L,I,G).

:- I\2=1, (T-P+360)\360>=(I+1)*90/G, orient(K,P), orient(L,T), not drop(N), same_opra(N,K,L,I,G).





% GUIDANCE FOR OPRA CONSTRAINTS

% Guide nondeterministic generation of orientation of objects

% when bounds are normal (no 360 modulo change)

:- I\2=1, J\2=1, K>L, (B-T+540)\360<=(((J-1)*90/G2)-((I+1)*90/G1)+360)\360, (((J-1)*90/G2)-((I+1)*90/G1)+360)\360<=(((J+1)*90/G2)-((I-1)*90/G1)+360)\360, orient(K,T), orient(L,B), single_diff_opra(M,K,L,I,G1), single_diff_opra(N,L,K,J,G2).

:- I\2=1, J\2=1, K>L, (B-T+540)\360>=(((J+1)*90/G2)-((I-1)*90/G1)+360)\360, (((J-1)*90/G2)-((I+1)*90/G1)+360)\360<=(((J+1)*90/G2)-((I-1)*90/G1)+360)\360, orient(K,T), orient(L,B), single_diff_opra(M,K,L,I,G1), single_diff_opra(N,L,K,J,G2).


% when bounds are reverse (there is 360 modulo change)

:- I\2=1, J\2=1, K>L, (B-T+540)\360<=(((J-1)*90/G2)-((I+1)*90/G1)+360)\360, (B-T+540)\360>=(((J+1)*90/G2)-((I-1)*90/G1)+360)\360, (((J-1)*90/G2)-((I+1)*90/G1)+360)\360>(((J+1)*90/G2)-((I-1)*90/G1)+360)\360, orient(K,T), orient(L,B), single_diff_opra(M,K,L,I,G1), single_diff_opra(N,L,K,J,G2).




% Guide nondeterministic generation of position of objects

:- I\2=1, Y<=0, ydist(K,L,Y), (360+T+(I-1)*90/G)\360<=180, (360+T+(I+1)*90/G)\360<=180, orient(L,T), single_diff_opra(N,K,L,I,G).
:- I\2=1, Y>=0, ydist(L,K,Y), (360+T+(I-1)*90/G)\360<=180, (360+T+(I+1)*90/G)\360<=180, orient(L,T), single_diff_opra(N,K,L,I,G).

:- I\2=1, Y>=0, ydist(K,L,Y), (360+T+(I-1)*90/G)\360>=180, (360+T+(I+1)*90/G)\360>=180, orient(L,T), single_diff_opra(N,K,L,I,G).
:- I\2=1, Y<=0, ydist(L,K,Y), (360+T+(I-1)*90/G)\360>=180, (360+T+(I+1)*90/G)\360>=180, orient(L,T), single_diff_opra(N,K,L,I,G).

:- I\2=1, X<=0, xdist(K,L,X), (270+T+(I-1)*90/G)\360>=180, (270+T+(I+1)*90/G)\360>=180, orient(L,T), single_diff_opra(N,K,L,I,G).
:- I\2=1, X>=0, xdist(L,K,X), (270+T+(I-1)*90/G)\360>=180, (270+T+(I+1)*90/G)\360>=180, orient(L,T), single_diff_opra(N,K,L,I,G).

:- I\2=1, X>=0, xdist(K,L,X), (90+T+(I-1)*90/G)\360>=180, (90+T+(I+1)*90/G)\360>=180, orient(L,T), single_diff_opra(N,K,L,I,G).
:- I\2=1, X<=0, xdist(L,K,X), (90+T+(I-1)*90/G)\360>=180, (90+T+(I+1)*90/G)\360>=180, orient(L,T), single_diff_opra(N,K,L,I,G).





#show xloc/2.
#show yloc/2.
#show orient/2.
