


% Infer unknown OPRA, position and distance relations

instantiate(P) :- point(K,P), toinfer_opra(K,L,G).
instantiate(R) :- point(L,R), toinfer_opra(K,L,G).

instantiate(P) :- point(K,P), toinfer_dist(K,L).
instantiate(R) :- point(L,R), toinfer_dist(K,L).

instantiate(P) :- point(K,P), toinfer_distrel(K,L).
instantiate(R) :- point(L,R), toinfer_distrel(K,L).

instantiate(P) :- point(K,P), toinfer_psi(K,L).
instantiate(R) :- point(L,R), toinfer_psi(K,L).

instantiate(P) :- point(K,P), toinfer_location(K).

generate_orient(L) :- toinfer_opra(K,L,G).

generate_orient(K) :- toinfer_orient(K).

generate_orient(L) :- toinfer_psi(K,L).



inferred_orient(K,A) :- orient(K,A), toinfer_orient(K).

inferred_location(K,X,Y) :- xloc(P,X), yloc(P,Y), point(K,P), toinfer_location(K).


compute_vector(P,R) :- point(K,P), point(L,R), toinfer_opra(K,L,G).

compute_vector(P,R) :- point(K,P), point(L,R), toinfer_psi(K,L).



same_loc(K,L) :- xdist(P,R,0), ydist(P,R,0), point(K,P), point(L,R), toinfer_psi(K,L).

inferred_delta(K,L,(B-T+360)\360) :- orient(K,B), orient(L,T), same_loc(K,L), toinfer_psi(K,L).

inferred_psi(K,L,(A-T+360)\360) :- orient(L,T), varphi(P,R,A), point(K,P), point(L,R), not same_loc(K,L), toinfer_psi(K,L).



compute_loc_diff(P,R) :- not compute_vector(P,R), point(K,P), point(L,R), toinfer_dist(K,L).

compute_loc_diff(P,R) :- not compute_vector(P,R), point(K,P), point(L,R), toinfer_distrel(K,L).

compute_dist(P,R) :- point(K,P), point(L,R), toinfer_dist(K,L).

compute_dist(P,R) :- point(K,P), point(L,R), toinfer_distrel(K,L).




% Infer OPRA relations


same_loc(K,L) :- xdist(P,R,0), ydist(P,R,0), point(K,P), point(L,R), toinfer_opra(K,L,G).


orient_diff(K,L,(B-T+360)\360) :- orient(K,B), orient(L,T), same_loc(K,L), toinfer_opra(K,L,G).

inferred_same_opra(K,L,Z/(90/G),G) :- Z\(180/G)=0, orient_diff(K,L,Z), same_loc(K,L), toinfer_opra(K,L,G). 

inferred_same_opra(K,L,1+2*(Z/(180/G)),G) :- Z\(180/G)!=0, orient_diff(K,L,Z), same_loc(K,L), toinfer_opra(K,L,G). 



angle_diff(K,L,(A-T+360)\360) :- orient(L,T), varphi(P,R,A), point(K,P), point(L,R), not same_loc(K,L), toinfer_opra(K,L,G).

inferred_diff_opra(K,L,D/(90/G),G) :- D\(180/G)=0, angle_diff(K,L,D), not same_loc(K,L), toinfer_opra(K,L,G).

inferred_diff_opra(K,L,1+2*(D/(180/G)),G) :- D\(180/G)!=0, angle_diff(K,L,D), not same_loc(K,L), toinfer_opra(K,L,G).



% Compute inferred distance relations

inferred_dist_sq(K,L,D) :- ndist(P,R,D), point(K,P), point(L,R), toinfer_dist(K,L).

inferred_distrel(K,L,U) :- D>=D1, D<=D2, lower_bound(U,D1), upper_bound(U,D2), ndist(P,R,D), point(K,P), point(L,R), toinfer_distrel(K,L).





#show inferred_orient/2.
#show inferred_location/3.
#show inferred_psi/3.
#show inferred_delta/3.
#show inferred_diff_opra/4.
#show inferred_same_opra/4.
#show inferred_dist_sq/3.
#show inferred_distrel/3.

