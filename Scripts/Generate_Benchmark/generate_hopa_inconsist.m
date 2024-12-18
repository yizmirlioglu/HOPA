function [outcome] = generate_hopa_inconsist()




    function outcome = writefile_asp(afolder, nobject, nconstraint, gr, carray, num_disjconst, num_disjx, ind_array, num_def, def_array, orgpair)

       if (num_disjconst > 0)   %(num_disjconst == 0 && num_def == 0)
          fnamex = cstrcat(afolder,"/n", num2str(nobject), "_c", num2str(nconstraint), "_grid", num2str(gr), "_disj", num2str(num_disjconst), "x", num2str(num_disjx), "_unsat.lp");
       elseif (num_def > 0)   %(num_disjconst == 0 && num_def == 0)
          fnamex = cstrcat(afolder,"/n", num2str(nobject), "_c", num2str(nconstraint), "_grid", num2str(gr), "_pres", num2str(num_def), "_unsat.lp");
       else
          fnamex = cstrcat(afolder,"/n", num2str(nobject), "_c", num2str(nconstraint), "_grid", num2str(gr), "_basic_unsat.lp");
       endif

       fptr = fopen(fnamex,"w");
       fprintf(fptr,"\n");
       fprintf(fptr,cstrcat("object(1..",num2str(nobject),").") );
       fprintf(fptr,"\n\n");

       for aobj=1:nobject
          fprintf(fptr,cstrcat("point(", num2str(aobj), ",", num2str(aobj), ").  ") );
       endfor    % for each object
       fprintf(fptr,"\n\n");

       %fprintf(fptr,cstrcat("xcoord(1..",num2str(gr),").\n") );
       %fprintf(fptr,cstrcat("ycoord(1..",num2str(gr),").\n\n") );
       %fprintf(fptr,"\n\n");


       if (num_def == 0)
          fprintf(fptr,cstrcat("constraint_id(1..", num2str(nconstraint), ").") );
       else
         for xindx=1:nconstraint
           if (ismember(xindx,def_array))
               fprintf(fptr,cstrcat("presumed_constraint_id(", num2str(xindx), "). ") );
           else
               fprintf(fptr,cstrcat("constraint_id(", num2str(xindx), "). ") );
           endif
         endfor
       endif

       fprintf(fptr,"\n\n");

       %fprintf(fptr,cstrcat("constraint_id(1..", num2str(nconstraint-num_def), ").") );
       %fprintf(fptr,"\n\n");
       %if (num_def > 0)
       %   fprintf(fptr,cstrcat("presumed_constraint_id(", num2str(nconstraint-num_def+1), "..", num2str(nconstraint), ").") );
       %   fprintf(fptr,"\n\n");
       %endif

        dconstr_ind = 1;
        for cindx=1:nconstraint
          aconstr = carray{cindx,1};
          if(ind_array(cindx) == 0)       % BASIC

             if (num_def > 0 && ismember(cindx,def_array) )   % cindx > nconstraint-num_def)     % presumed constraint
                fprintf(fptr,"presumed_");
             endif

             if (aconstr(1) == 1)     % same opra
                fprintf(fptr,cstrcat("same_opra(", num2str(cindx), ",", num2str(aconstr(3)), ",", num2str(aconstr(2)), ",", num2str(aconstr(4)), ",", num2str(aconstr(5)), "). \n" ) );
             elseif (aconstr(1) == 2)     % one-sided diff opra
                   fprintf(fptr,cstrcat("one_sided_diff_opra(", num2str(cindx), ",", num2str(aconstr(2)), ",", num2str(aconstr(3)), ",", num2str(aconstr(4)), ",", num2str(aconstr(5)),  "). \n" ) );
             elseif (aconstr(1) == 3)     % diff opra
                   fprintf(fptr,cstrcat("diff_opra(", num2str(cindx), ",", num2str(aconstr(2)), ",", num2str(aconstr(3)), ",", num2str(aconstr(4)), ",", num2str(aconstr(5)), ",", num2str(aconstr(6)), "). \n" ) );
             elseif (aconstr(1) == 4)     % qualitative dist
                   fprintf(fptr,cstrcat("dist_rel(", num2str(cindx), ",", num2str(aconstr(2)), ",", num2str(aconstr(3)), ",", dist_name{aconstr(4)}, "). \n" ) );
             elseif (aconstr(1) == 5)     % 4:orient exact
                   fprintf(fptr,cstrcat("direction(", num2str(cindx), ",", num2str(aconstr(2)), ",", num2str(aconstr(3)), "). \n" ) );
             elseif (aconstr(1) == 6)     % 5:orient range
                   fprintf(fptr,cstrcat("direction_range(", num2str(cindx), ",", num2str(aconstr(2)), ",", num2str(aconstr(3)), ",", num2str(aconstr(4)), "). \n" ) );
             elseif (aconstr(1) == 7)     % orient vector
                   fprintf(fptr,cstrcat("direction(", num2str(cindx), ",", num2str(aconstr(2)), ",", num2str(aconstr(3)), ",", num2str(aconstr(4)), "). \n" ) );
             elseif (aconstr(1) == 8)     % psi exact
                   fprintf(fptr,cstrcat("psi_precise(", num2str(cindx), ",", num2str(aconstr(2)), ",", num2str(aconstr(3)), ",", num2str(aconstr(4)), "). \n" ) );
             elseif (aconstr(1) == 9)     % psi range
                   fprintf(fptr,cstrcat("psi_range(", num2str(cindx), ",", num2str(aconstr(2)), ",", num2str(aconstr(3)), ",", num2str(aconstr(4)), ",", num2str(aconstr(5)), "). \n" ) );
             elseif (aconstr(1) == 10)     % delta exact
                   fprintf(fptr,cstrcat("delta_precise(", num2str(cindx), ",", num2str(aconstr(3)), ",", num2str(aconstr(2)), ",", num2str(aconstr(4)), "). \n" ) );
             elseif (aconstr(1) == 11)     % delta range
                   fprintf(fptr,cstrcat("delta_range(", num2str(cindx), ",", num2str(aconstr(3)), ",", num2str(aconstr(2)), ",", num2str(aconstr(4)), ",", num2str(aconstr(5)), "). \n" ) );
             elseif (aconstr(1) == 12)     % exact location
                   fprintf(fptr,cstrcat("location(", num2str(cindx), ",", num2str(aconstr(2)), ",", num2str(aconstr(3)), ",", num2str(aconstr(4)), "). \n" ) );
             elseif (aconstr(1) == 13)     % distance range
                   fprintf(fptr,cstrcat("distance_range(", num2str(cindx), ",", num2str(aconstr(2)), ",", num2str(aconstr(3)), ",", num2str(aconstr(4)), ",", num2str(aconstr(5)), "). \n" ) );
             elseif (aconstr(1) == 14)     % exact distance
                   fprintf(fptr,cstrcat("distance(", num2str(cindx), ",", num2str(aconstr(2)), ",", num2str(aconstr(3)), ",", num2str(aconstr(4)), "). \n" ) );

             endif
                %fprintf(fptr,cstrcat("dist_rel(", num2str(aconstr(1)), ",", num2str(aconstr(2)), ",", dist_name{aconstr(7)}, "). \n" ) );

          else              % disjunctive

             fprintf(fptr,cstrcat("disj_index(", num2str(cindx), ",1..", num2str(num_disjx), "). \n" ) );
             fprintf(fptr,cstrcat("original_pair(", num2str(cindx), ",", num2str(orgpair{dconstr_ind}(1)), ",", num2str(orgpair{dconstr_ind}(2)), "). \n" ) );
             for (dindx=1:num_disjx)
                dconstr = aconstr{dindx};

               if (dconstr(1) == 1)     % same opra
                   fprintf(fptr,cstrcat("disjrelation(", num2str(cindx), ",", num2str(dindx), ",",  num2str(dconstr(1)), ",", num2str(dconstr(3)), ",", num2str(dconstr(2)), ",", num2str(dconstr(4)), ",", num2str(dconstr(5)), "). \n" ) );
               elseif (dconstr(1) == 2)     % one-sided diff opra
                   fprintf(fptr,cstrcat("disjrelation(", num2str(cindx), ",", num2str(dindx), ",",  num2str(dconstr(1)), ",", num2str(dconstr(2)), ",", num2str(dconstr(3)), ",", num2str(dconstr(4)), ",", num2str(dconstr(5)), "). \n" ) );
               elseif (dconstr(1) == 3)     % diff opra
                   fprintf(fptr,cstrcat("disjrelation(", num2str(cindx), ",", num2str(dindx), ",",  num2str(dconstr(1)), ",", num2str(dconstr(2)), ",", num2str(dconstr(3)), ",", num2str(dconstr(4)), ",", num2str(dconstr(5)), ",", num2str(dconstr(6)), "). \n" ) );
               elseif (dconstr(1) == 4)     % qualitative dist
                   fprintf(fptr,cstrcat("disjrelation(", num2str(cindx), ",", num2str(dindx), ",",  num2str(dconstr(1)), ",", num2str(dconstr(2)), ",", num2str(dconstr(3)), ",", dist_name{dconstr(4)}, "). \n" ) );
               elseif (dconstr(1) == 5)     % 4:orient exact
                   fprintf(fptr,cstrcat("disjrelation(", num2str(cindx), ",", num2str(dindx), ",",  num2str(dconstr(1)), ",", num2str(dconstr(2)), ",", num2str(dconstr(3)), "). \n" ) );
               elseif (dconstr(1) == 6)     % 5:orient range
                   fprintf(fptr,cstrcat("disjrelation(", num2str(cindx), ",", num2str(dindx), ",",  num2str(dconstr(1)), ",", num2str(dconstr(2)), ",", num2str(dconstr(3)), ",", num2str(dconstr(4)), "). \n" ) );
               elseif (dconstr(1) == 7)     % orient vector
                   fprintf(fptr,cstrcat("disjrelation(", num2str(cindx), ",", num2str(dindx), ",",  num2str(dconstr(1)), ",", num2str(dconstr(2)), ",", num2str(dconstr(3)), ",", num2str(dconstr(4)), "). \n" ) );
               elseif (dconstr(1) == 8)     % psi exact
                   fprintf(fptr,cstrcat("disjrelation(", num2str(cindx), ",", num2str(dindx), ",",  num2str(dconstr(1)), ",", num2str(dconstr(2)), ",", num2str(dconstr(3)), ",", num2str(dconstr(4)), "). \n" ) );
               elseif (dconstr(1) == 9)     % psi range
                   fprintf(fptr,cstrcat("disjrelation(", num2str(cindx), ",", num2str(dindx), ",",  num2str(dconstr(1)), ",", num2str(dconstr(2)), ",", num2str(dconstr(3)), ",", num2str(dconstr(4)), ",", num2str(dconstr(5)), "). \n" ) );
               elseif (dconstr(1) == 10)     % delta exact
                   fprintf(fptr,cstrcat("disjrelation(", num2str(cindx), ",", num2str(dindx), ",",  num2str(dconstr(1)), ",", num2str(dconstr(3)), ",", num2str(dconstr(2)), ",", num2str(dconstr(4)), "). \n" ) );
               elseif (dconstr(1) == 11)     % delta range
                   fprintf(fptr,cstrcat("disjrelation(", num2str(cindx), ",", num2str(dindx), ",",  num2str(dconstr(1)), ",", num2str(dconstr(3)), ",", num2str(dconstr(2)), ",", num2str(dconstr(4)), ",", num2str(dconstr(5)), "). \n" ) );
               elseif (dconstr(1) == 12)     % exact location
                   fprintf(fptr,cstrcat("disjrelation(", num2str(cindx), ",", num2str(dindx), ",",  num2str(dconstr(1)), ",", num2str(dconstr(2)), ",", num2str(dconstr(3)), ",", num2str(dconstr(4)), "). \n" ) );
               elseif (dconstr(1) == 13)     % distance range
                   fprintf(fptr,cstrcat("disjrelation(", num2str(cindx), ",", num2str(dindx), ",",  num2str(dconstr(1)), ",", num2str(dconstr(2)), ",", num2str(dconstr(3)), ",", num2str(dconstr(4)), ",", num2str(dconstr(5)), "). \n" ) );
               elseif (dconstr(1) == 14)     % exact distance
                   fprintf(fptr,cstrcat("disjrelation(", num2str(cindx), ",", num2str(dindx), ",",  num2str(dconstr(1)), ",", num2str(dconstr(2)), ",", num2str(dconstr(3)), ",", num2str(dconstr(4)), "). \n" ) );

                endif
  
             endfor
             dconstr_ind = dconstr_ind+1;

          endif  % basic or disjunctive
       endfor

       fclose(fptr);
       outcome = 1;
       
    endfunction




  dist_name = {"on", "verynear", "near", "medium", "far", "veryfar"};
  
  startindex = 1;
  num_samples = 100;

    % Number of objects, constraints, grid size
  num_objects = [40, 60, 80, 100, 120];  %,[10, 20, 40, 60]
  num_constr = {[40,60], [60, 100], [80, 120], [100,150], [120,150]};
  num_grid = { {[50,100],[50,100]}, {[50,100],[50,100]}, {[50,100],[50]}, {[50],[50]}, {[50],[50]} };

  disj_objects = [40, 60, 80];   %, 100, 200];  %, 50, 100, 200, 400, 500];
  disj_constr = {[40, 60], [60], [80]};
  disj_grid = {[50], [50], [50]}; 

  num_disjunc = [2,4,8];
  disjunc_ratio =  [0.2, 0.4];  % [0.1, 0.2];
  max_disjunc = 8;

  def_objects = [10, 20, 40];   %, 100, 200];  %, 50, 100, 200, 400, 500];
  def_constr = {[10, 20], [20], [40]};
  def_grid = {[50,100], [50], [50]}; 

  def_ratio =  [0.1, 0.2];  % [0.1, 0.2];

  dist_range_min = 0.1;  % param for max distance range type rel
  dist_range_max = 0.4;  % param for max distance range type rel
  orient_min_range = 15;  % param for psi or delta range type rel
  orient_max_range = 45;  % param for psi or delta range type rel


  % Distance relations 
  distrel = 1:6;    % 0:5
  dist_lower_ratio = [-0.1, 0, 0.10, 0.20, 0.40, 0.70];     % exclusive lower bounds
  dist_upper_ratio = [0, 0.10, 0.20, 0.40, 0.70, 1.00];     % inclusive upper bounds

  len_dist = length(distrel);

  
  aspfolder = "./ASP";

  consist_folder = "/Consistent";
  inconsist_folder = "/Inconsistent";

  sample_folder = "/Sample";

  basic_folder = "/Basic";
  disj_folder = "/Disjunctive";
  default_folder = "/Presumed";

  sleep_int = 2.0;
  sleep_sample = 10.0;
  
  len_objects = length(num_objects);
  len_disj_objects = length(disj_objects);
  len_disjunc = length(num_disjunc);
  len_ratio = length(disjunc_ratio);

   

  printf("------ RANDOM RELATION SAMPLES----- \n");
  printf("Finished Sample ");
  
  for m=1:num_samples
     sample_asp_folder = cstrcat(aspfolder,inconsist_folder,sample_folder,num2str(m+startindex-1));

     res1a = mkdir(sample_asp_folder);
     pause(sleep_int);

     sample_basp_folder = cstrcat(sample_asp_folder,basic_folder);

     sample_dasp_folder = cstrcat(sample_asp_folder,disj_folder);

     res3a = mkdir(sample_basp_folder);
     pause(sleep_int);
     res3c = mkdir(sample_dasp_folder);
     pause(sleep_int);
     
     for n=1:len_objects
        nobj = num_objects(n); 
        len_constr = length(num_constr{n});

        for k=1:len_constr
           nconstr = num_constr{n}(k);  
           constr_grid = (num_grid{n}){k};
           len_grid = length(constr_grid);   
      
           % check whether is in disjunctive constraints
           is_in_disj = 0;
           nindex = find(disj_objects == nobj);
           if (length(nindex) > 0)
               constr_ind = find(disj_constr{nindex} == nconstr);
               if (length(constr_ind) > 0)  
                     is_in_disj = 1;
               endif
           endif


          for z=1:len_grid

            gridx = constr_grid(z);
            ascale = gridx * sqrt(2);
            dist_lower_bounds = round(dist_lower_ratio * ascale);     % exclusive lower bounds
            dist_upper_bounds = round(dist_upper_ratio * ascale);     % exclusive lower bounds

            % check whether is in presumed constraints
            is_in_def = 0;
            defindex = find(def_objects == nobj);
            if (length(defindex) > 0)
                defconstr_ind = find(def_constr{defindex} == nconstr);
                if  ismember(nconstr, def_constr{defindex})  &&  ismember(gridx, def_grid{defindex})
                     is_in_def = 1;
                endif
            endif

            constr_type = randi(14,1,nconstr);  
            % 1: same opra, 2: one-sided diff opra 3: diff opra 4: qualit dist, 5:orient exact, 6:orient range, 7:orient vector, 8.psi exact, 9:psi range, 10:delta exact, 11. delt range, 12:exact location  13: distance range, 14: exact distance

             % constraint pairs
            trg_obj = randi(nobj,1,nconstr); 
            ref_obj = randi(nobj,1,nconstr); 

            constr_array = cell(nconstr,1);   

            for t=1:nconstr

               targt = trg_obj(1,t);
               refern = ref_obj(1,t);

               granul = 1 + randi(5);  

               if (targt == refern) && (refern > 1)    % target and ref objects cannot be the same
                  refern = refern - 1;
                  ref_obj(1,t) = ref_obj(1,t) + 1;
               elseif (targt == refern) && (refern == 1)
                  targt = 2;
                  trg_obj(1,t) = 2;
               endif

               if (constr_type(1,t) == 1)   % same opra
                  constr_array{t,1} = [ 1, targt, refern, randi(4*granul)-1, granul ];
               elseif (constr_type(1,t) == 2) % one sided diff opra
                  constr_array{t,1} = [ 2, targt, refern, randi(4*granul)-1, granul ];
               elseif (constr_type(1,t) == 3) % diff opra
                  constr_array{t,1} = [ 3, targt, refern, randi(4*granul)-1, randi(4*granul)-1, granul ];
               elseif (constr_type(1,t) == 4) % : qualit dist
                  constr_array{t,1} = [ 4, targt, refern, 1+randi(len_dist-1) ];
               elseif (constr_type(1,t) == 5) % :orient exact
                  constr_array{t,1} = [ 5, targt, randi(360)-1 ];
               elseif (constr_type(1,t) == 6) % :orient range
                  randangle = randi(360-orient_max_range,1,2) - 1;
                  arangex = orient_min_range + randi(orient_max_range-orient_min_range);
                  constr_array{t,1} = [ 6, targt, max(randangle-arangex,0), min(randangle+arangex,359) ];
               elseif (constr_type(1,t) == 7) % orient vector
                  arangex = randi(nobj,1,2);
                  if (arangex(1) == arangex(2)) && (arangex(1) > 1)    % target and ref objects cannot be the same
                     arangex(1) = arangex(1) - 1;
                  elseif (arangex(1) == arangex(2)) && (arangex(1) == 1)
                     arangex(1) = 2;
                  endif
                  constr_array{t,1} = [ 7, targt, min(arangex), max(arangex) ];
               elseif (constr_type(1,t) == 8) % psi exact
                  constr_array{t,1} = [ 8, targt, refern, randi(360)-1 ];
               elseif (constr_type(1,t) == 9) % psi range
                  randangle = randi(360-orient_max_range,1,2) - 1;
                  arangex = orient_min_range + randi(orient_max_range-orient_min_range);
                  constr_array{t,1} = [9, targt, refern, max(randangle-arangex,0), min(randangle+arangex,359) ];
               elseif (constr_type(1,t) == 10) % delta exact
                  constr_array{t,1} = [ 10, targt, refern, randi(360)-1 ];
               elseif (constr_type(1,t) == 11) % delta range
                  randangle = randi(360-orient_max_range,1,2) - 1;
                  arangex = orient_min_range + randi(orient_max_range-orient_min_range);
                  constr_array{t,1} = [11, targt, refern, max(randangle-arangex,0), min(randangle+arangex,359) ];
               elseif (constr_type(1,t) == 12) % exact location
                  constr_array{t,1} = [ 12, targt, randi(gridx), randi(gridx) ];
               elseif (constr_type(1,t) == 13) % distance range
                  arangex = randi(2*gridx*gridx,1,2);
                  constr_array{t,1} = [ 13, targt, refern, min(arangex), max(arangex) ];
               elseif (constr_type(1,t) == 14) % exact distance
                  constr_array{t,1} = [ 14, targt, refern, randi(2*gridx*gridx) ];
               endif

             endfor  % number of basic constraints

              res6 = writefile_asp(sample_basp_folder,nobj,nconstr,gridx,constr_array,0,0,zeros(nconstr,1),0,[],[]);
              pause(sleep_int);


              if (is_in_disj == 1) 
                if ismember(gridx, disj_grid{nindex})
 
                 for rind=1:len_ratio
                   ndisj_constr = round(disjunc_ratio(rind)*nconstr);

                  constr_ind = zeros(nconstr,1);   % 0 for basic
                  disj_constr_ind = sort(randperm(nconstr,ndisj_constr));   % choose random constraints to convert to disjunctive
                  constr_ind(disj_constr_ind) = 1;   % choose random constraints to convert to disjunctive

                  for disj_ind=1:len_disjunc
                     ndisjunc = num_disjunc(disj_ind);

                     target_objs = randi(nobj,ndisj_constr,ndisjunc-1); 
                     refer_objs = randi(nobj,ndisj_constr,ndisjunc-1); 
                     dconstr_type = randi(14,ndisj_constr,ndisjunc-1);  
                     % 1: same opra, 2: one-sided diff opra 3: diff opra 4: qualit dist, 5:orient exact, 6:orient range, 7:orient vector, 8.psi exact, 9:psi range, 10:delta exact, 11. delt range, 12:exact location, 13: distance range, 14: exact distance

                     disj_array = constr_array;   %cell(ndisjunc-1,nsoft);
                     orig_pair = {};  

                     for cind=1:ndisj_constr
                        constr_id = disj_constr_ind(cind);
                        disj_array{constr_id,1} = cell(ndisjunc,1);
                        disj_array{constr_id,1}{1} = constr_array{constr_id,1};

                        orig_pair{cind} = [trg_obj(1,constr_id), ref_obj(1,constr_id)]; 

                        for dindex=1:(ndisjunc-1)

                          trgd = target_objs(cind,dindex);
                          refd = refer_objs(cind,dindex);
                          granul = 1 + randi(5);  

                          if (trgd == refd) && (refd > 1)    % target and ref objects cannot be the same
                              refd = refd - 1;
                          elseif (trgd == refd) && (refd == 1)
                              trgd = 2;
                          endif

                          if (dconstr_type(cind,dindex) == 1)   % same opra
                             disj_array{constr_id,1}{dindex+1} = [ 1, trgd, refd, randi(4*granul)-1, granul ];
                          elseif (dconstr_type(cind,dindex) == 2) % one sided diff opra
                             disj_array{constr_id,1}{dindex+1} = [ 2, trgd, refd, randi(4*granul)-1, granul ];
                          elseif (dconstr_type(cind,dindex) == 3) % diff opra
                             disj_array{constr_id,1}{dindex+1} = [ 3, trgd, refd, randi(4*granul)-1, randi(4*granul)-1, granul ];
                          elseif (dconstr_type(cind,dindex) == 4) % : qualit dist
                             disj_array{constr_id,1}{dindex+1} = [ 4, trgd, refd, 1+randi(len_dist-1) ];
                          elseif (dconstr_type(cind,dindex) == 5) % :orient exact
                             disj_array{constr_id,1}{dindex+1} = [ 5, trgd, randi(360)-1 ];
                          elseif (dconstr_type(cind,dindex) == 6) % :orient range
                              randangle = randi(360-orient_max_range,1,2) - 1;
                              arangex = orient_min_range + randi(orient_max_range-orient_min_range);
                              disj_array{constr_id,1}{dindex+1} = [ 6, trgd, max(randangle-arangex,0), min(randangle+arangex,359) ];
                          elseif (dconstr_type(cind,dindex) == 7) % orient vector
                             arangex = randi(nobj,1,2);
                             if (arangex(1) == arangex(2)) && (arangex(1) > 1)    % target and ref objects cannot be the same
                               arangex(1) = arangex(1) - 1;
                             elseif (arangex(1) == arangex(2)) && (arangex(1) == 1)
                               arangex(1) = 2;
                             endif
                              disj_array{constr_id,1}{dindex+1} = [ 7, trgd, min(arangex), max(arangex) ];
                          elseif (dconstr_type(cind,dindex) == 8) % psi exact
                              disj_array{constr_id,1}{dindex+1} = [ 8, trgd, refd, randi(360)-1 ];
                          elseif (dconstr_type(cind,dindex) == 9) % psi range
                              randangle = randi(360-orient_max_range,1,2) - 1;
                              arangex = orient_min_range + randi(orient_max_range-orient_min_range);
                              disj_array{constr_id,1}{dindex+1} = [ 9, trgd, refd, max(randangle-arangex,0), min(randangle+arangex,359) ];
                          elseif (dconstr_type(cind,dindex) == 10) % delta exact
                             disj_array{constr_id,1}{dindex+1} = [ 10, trgd, refd, randi(360)-1 ];
                          elseif (dconstr_type(cind,dindex) == 11) % delta range
                              randangle = randi(360-orient_max_range,1,2) - 1;
                              arangex = orient_min_range + randi(orient_max_range-orient_min_range);
                              disj_array{constr_id,1}{dindex+1} = [ 11, trgd, refd, max(randangle-arangex,0), min(randangle+arangex,359) ];
                          elseif (dconstr_type(cind,dindex) == 12) % exact location
                             disj_array{constr_id,1}{dindex+1} = [ 12, trgd, randi(gridx), randi(gridx) ];
                          elseif (dconstr_type(cind,dindex) == 13) % distance range
                             arangex = randi(2*gridx*gridx,1,2);
                             disj_array{constr_id,1}{dindex+1} = [ 13, trgd, refd, min(arangex), max(arangex) ];
                          elseif (dconstr_type(cind,dindex) == 14) % exact distance
                             disj_array{constr_id,1}{dindex+1} = [ 14, trgd, refd, randi(2*gridx*gridx) ];
                          endif

                        endfor    % for each disjunct
                     endfor    % for each disjunctive constraint


                     res9 = writefile_asp(sample_dasp_folder,nobj,nconstr,gridx,disj_array,ndisj_constr,ndisjunc,constr_ind,0,[],orig_pair); %{rind,disj_ind}
                      pause(sleep_int);

                   endfor    % num disjuncs
                 endfor   %  disjunctive ratio

              endif  % if in disj grid
             endif  % if in disjunc size

          endfor    % grids

       endfor   %  num constraints

     endfor       % num objects

     printf( cstrcat(num2str(m+startindex-1), " ") );            %printf(cstrcat("Finished Sample ", num2str(m+startindex-1)) );
     pause(sleep_sample);

   endfor       % samples 

  
  printf("------ Produced Inconsistent Instances----- \n");
  outcome =1;

  
endfunction
