

% Presumed information about position and orientation of objects

{violated(N)}1 :- presumed_constraint_id(N).


%  Optimization:  minimize the number of violated presumed constraints

:~ violated(N).  [1@2,N]



% If presumed HOPA constraint is not violated, it becomes a hard HOPA constraint

diff_opra(N,K,L,I,J,G) :- not violated(N), presumed_diff_opra(N,K,L,I,J,G).

same_opra(N,K,L,I,G) :- not violated(N), presumed_same_opra(N,K,L,I,G).

one_sided_diff_opra(N,K,L,I,G) :- not violated(N), presumed_one_sided_diff_opra(N,K,L,I,G).

location(N,K,X,Y) :- not violated(N), presumed_location(N,K,X,Y).

direction(N,K,A) :- not violated(N), presumed_direction(N,K,A).

direction(N,K,P,R) :- not violated(N), presumed_direction(N,K,P,R).

direction_range(N,K,A1,A2) :- not violated(N), presumed_direction_range(N,K,A1,A2).

psi_precise(N,K,L,A) :- not violated(N), presumed_psi_precise(N,K,L,A).

psi_range(N,K,L,A1,A2) :- not violated(N), presumed_psi_range(N,K,L,A1,A2).

dist_rel(N,K,L,U) :- not violated(N), presumed_dist_rel(N,K,L,U).

distance(N,K,L,D) :- not violated(N), presumed_distance(N,K,L,D).

distance_range(N,K,L,D1,D2) :- not violated(N), presumed_distance_range(N,K,L,D1,D2).




% Presumed Disjunctive constraints

% 1: same opra, 2: one-sided diff opra 3: diff opra 4: qualit dist, 5:orient exact, 6:orient range, 7:orient vector, 8.psi exact, 9:psi range, 10:delta exact, 11. delt range, 12:exact location  13: distance range, 14: exact distance


exist_pdisj(N) :- not violated(N), presumed_disj_index(N,D). 


1{pchosen(N,D): presumed_disj_index(N,D)}1 :- not violated(N), exist_pdisj(N).


same_opra(N,K,L,I,G) :- not violated(N), pchosen(N,D), presumed_disjrelation(N,D,1,K,L,I,G).

one_sided_diff_opra(N,K,L,I,G) :- not violated(N), pchosen(N,D), presumed_disjrelation(N,D,2,K,L,I,G).

diff_opra(N,K,L,I,J,G) :- not violated(N), pchosen(N,D), presumed_disjrelation(N,D,3,K,L,I,J,G).

dist_rel(N,K,L,U) :- not violated(N), pchosen(N,D), presumed_disjrelation(N,D,4,K,L,U).

direction(N,K,A) :- not violated(N), pchosen(N,D), presumed_disjrelation(N,D,5,K,A).

direction_range(N,K,A1,A2) :- not violated(N), pchosen(N,D), presumed_disjrelation(N,D,6,K,A1,A2).

direction(N,K,P,R) :- not violated(N), pchosen(N,D), presumed_disjrelation(N,D,7,K,P,R).

psi_precise(N,K,L,A) :- not violated(N), pchosen(N,D), presumed_disjrelation(N,D,8,K,L,A).

psi_range(N,K,L,A1,A2) :- not violated(N), pchosen(N,D), presumed_disjrelation(N,D,9,K,L,A1,A2).

delta_precise(N,K,L,A) :- not violated(N), pchosen(N,D), presumed_disjrelation(N,D,10,K,L,A).

delta_range(N,K,L,A1,A2) :- not violated(N), pchosen(N,D), presumed_disjrelation(N,D,11,K,L,A1,A2).

location(N,K,X,Y) :- not violated(N), pchosen(N,D), presumed_disjrelation(N,D,12,K,X,Y).

distance_range(N,K,L,U1,U2) :- not violated(N), pchosen(N,D), presumed_disjrelation(N,D,13,K,L,U1,U2).

distance(N,K,L,U) :- not violated(N), pchosen(N,D), presumed_disjrelation(N,D,14,K,L,U).



#show violated/1.
#show pchosen/2.
