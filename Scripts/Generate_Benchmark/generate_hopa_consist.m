function [outcome] = generate_hopa_consist()
             
    function [is_same, sector1, sector2, angle_out1, angle_out2] = find_opra(trg_x, trg_y, trg_dir, ref_x, ref_y, ref_dir, gr)

    quadr = 180/gr;

    if ((trg_x != ref_x) || (trg_y != ref_y))  

       is_same = 0;

       xvec = trg_x-ref_x;
       yvec = trg_y-ref_y;
       if (xvec != 0)
          angle_deg = atand(abs(yvec)/abs(xvec));
       else
          angle_deg = 90;
       endif

       if ((yvec >= 0) && (xvec >= 0)) 
          angle_vec = angle_deg;
       elseif ((yvec >= 0) && (xvec < 0)) 
          angle_vec = 180 - angle_deg;
       elseif ((yvec < 0) && (xvec >= 0))  
          angle_vec = 360 - angle_deg;
       else
          angle_vec = 180 + angle_deg;
       endif 

       if (angle_vec < 180) 
          angle_rev = angle_vec + 180;
       else
          angle_rev = angle_vec - 180;
       endif

       angle_diff1 = angle_vec - ref_dir;
       if (angle_diff1 < 0) 
          angle_diff1 = angle_diff1 + 360;
       endif

       sect_coef1 = floor(angle_diff1/quadr);
       sect_mod1 = mod(angle_diff1, quadr);
      
       if (sect_mod1 > 0) 
          sector1 = 2*sect_coef1 + 1;
       else
          sector1 = 2*sect_coef1;
       endif

       angle_diff2 = angle_rev - trg_dir;
       if (angle_diff2 < 0) 
          angle_diff2 = angle_diff2 + 360;
       endif

       sect_coef2 = floor(angle_diff2/quadr);
       sect_mod2 = mod(angle_diff2, quadr);
      
       if (sect_mod2 > 0) 
          sector2 = 2*sect_coef2 + 1;
       else
          sector2 = 2*sect_coef2;
       endif

    else  % same opra

       is_same = 1;

       angle_diff1 = trg_dir - ref_dir;
       if (angle_diff1 < 0) 
          angle_diff1 = angle_diff1 + 360;
       endif

       sect_coef1 = floor(angle_diff1/quadr);
       sect_mod1 = mod(angle_diff1, quadr);
      
       if (sect_mod1 > 0) 
          sector1 = 2*sect_coef1 + 1;
       else
          sector1 = 2*sect_coef1;
       endif

       sector2 =  mod(4*gr-sector1, 4*gr);

       angle_diff2 = ref_dir - trg_dir;
       if (angle_diff2 < 0) 
          angle_diff2 = angle_diff2 + 360;
       endif

       %sect_coef2 = floor(angle_diff2/quadr);
       %sect_mod2 = mod(angle_diff2, quadr);
      
       %if (sect_mod2 > 0) 
       %   sector2 = 2*sect_coef2 + 1;
       %else
       %   sector2 = 2*sect_coef2;
       %endif

    endif

    angle_out1 = round(angle_diff1);
    angle_out2 = round(angle_diff2);

  
    endfunction





   function outcome = writefile_asp(afolder, nobject, nconstraint, gr, carray, num_disjconst, num_disjx, ind_array, num_def, def_array, numexp, orgpair)

       if (num_disjconst > 0)   %(num_disjconst == 0 && num_def == 0)
          fnamex = cstrcat(afolder,"/n", num2str(nobject), "_c", num2str(nconstraint), "_grid", num2str(gr), "_disj", num2str(num_disjconst), "x", num2str(num_disjx), "_sat.lp");
       elseif (num_def > 0)   %(num_disjconst == 0 && num_def == 0)
          fnamex = cstrcat(afolder,"/n", num2str(nobject), "_c", num2str(nconstraint), "_grid", num2str(gr), "_pres", num2str(num_def), "_sat.lp");
       elseif (numexp > 0)   %(num_disjconst == 0 && num_def == 0)
          fnamex = cstrcat(afolder,"/n", num2str(nobject), "_c", num2str(nconstraint), "_grid", num2str(gr), "_exp", num2str(numexp), "_sat.lp");
       else
          fnamex = cstrcat(afolder,"/n", num2str(nobject), "_c", num2str(nconstraint), "_grid", num2str(gr), "_basic_sat.lp");
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
             elseif (aconstr(1) == 5)     % orient exact
                   fprintf(fptr,cstrcat("direction(", num2str(cindx), ",", num2str(aconstr(2)), ",", num2str(aconstr(3)), "). \n" ) );
             elseif (aconstr(1) == 6)     % orient range
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


  def_objects = [40, 60, 80];   %, 100, 200];  %, 50, 100, 200, 400, 500];
  def_constr = {[40, 60], [60], [80]};
  def_grid = {[50], [50], [50]}; 

  def_ratio =  [0.1, 0.2];  % [0.1, 0.2];

  explain_objects = [40, 60];   %, 100, 200];  %, 50, 100, 200, 400, 500];
  explain_constr = {[40], [60]};
  explain_grid = {[50], [50], [50]}; 

  explain_ratio =  [0.1, 0.2];  % [0.1, 0.2];

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
  explain_folder = "/Explain";

  sleep_int = 2.0;
  sleep_sample = 12.0;
  
  len_objects = length(num_objects);
  len_disj_objects = length(disj_objects);
  len_disjunc = length(num_disjunc);
  len_ratio = length(disjunc_ratio);

   

  printf("------ RANDOM OBJECT SAMPLES----- \n");
  printf("Finished Sample ");
  
  for m=1:num_samples
     sample_asp_folder = cstrcat(aspfolder,consist_folder,sample_folder,num2str(m+startindex-1));

     res1a = mkdir(sample_asp_folder);
     pause(sleep_int);

     sample_basp_folder = cstrcat(sample_asp_folder,basic_folder);
     sample_dasp_folder = cstrcat(sample_asp_folder,disj_folder);
     sample_fasp_folder = cstrcat(sample_asp_folder,default_folder);
     sample_easp_folder = cstrcat(sample_asp_folder,explain_folder);

     res3a = mkdir(sample_basp_folder);
     pause(sleep_int);
     res3b = mkdir(sample_fasp_folder);
     pause(sleep_int);
     res3c = mkdir(sample_dasp_folder);
     pause(sleep_int);
     res5c = mkdir(sample_easp_folder);
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
               %gr_ind = find(disj_granul{nindex} == granul);
               if (length(constr_ind) > 0)    %  &&  length(gr_ind) > 0)
                     is_in_disj = 1;
               endif
           endif

        for z=1:len_grid

            gridx = constr_grid(z);
            ascale = gridx * sqrt(2);
            dist_lower_bounds = round(dist_lower_ratio * ascale);     % exclusive lower bounds
            dist_upper_bounds = round(dist_upper_ratio * ascale);     % exclusive lower bounds

            dist_range1 = round(dist_range_min * 2 * gridx * gridx);  %  
            dist_range2 = round(dist_range_max * 2 * gridx * gridx);  %  

            % check whether is in presumed constraints
            is_in_def = 0;
            defindex = find(def_objects == nobj);
            if (length(defindex) > 0)
                defconstr_ind = find(def_constr{defindex} == nconstr);
                if  ismember(nconstr, def_constr{defindex})  &&  ismember(gridx, def_grid{defindex})
                     is_in_def = 1;
                endif
            endif

            % check whether is in explain constraints
            is_in_explain = 0;
            if (ismember(nobj, explain_objects))   
                explainindex = find(explain_objects == nobj);
                if  ismember(nconstr, explain_constr{explainindex})  &&  ismember(gridx, explain_grid{explainindex})
                     is_in_explain = 1;
                endif
            endif


           % create objects
           obj_locx = randi(gridx,1,nobj); 
           obj_locy = randi(gridx,1,nobj);
           obj_dir = randi(360,1,nobj) - 1;  
          

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

               [issame, sect1, sect2, ang1, ang2] = find_opra(obj_locx(targt), obj_locy(targt), obj_dir(targt), obj_locx(refern), obj_locy(refern), obj_dir(refern), granul);

               if (constr_type(1,t) < 4)   % opra type rel
                  if (issame == 1)   % same opra
                     constr_array{t,1} = [ 1, targt, refern, sect1, granul ];
                  elseif (constr_type(1,t) == 1) % same opra  but not same location
                     constr_array{t,1} = [ 2, targt, refern, sect1, granul ];
                  elseif (constr_type(1,t) == 2) % one sided diff opra
                     constr_array{t,1} = [ 2, targt, refern, sect1, granul ];
                  elseif (constr_type(1,t) == 3) % diff opra
                     constr_array{t,1} = [ 3, targt, refern, sect1, sect2, granul ];
                  endif
               elseif (constr_type(1,t) == 4) % : qualit dist
                    % find distance relation
                  numdist = sqrt( (obj_locx(targt)-obj_locx(refern))^2 + (obj_locy(targt)-obj_locy(refern))^2 );
                  ind2 = numdist <= dist_upper_bounds; 
                  distind = find(ind2);   
                  constr_array{t,1} = [ 4, targt, refern, distind(1) ];
               elseif (constr_type(1,t) == 5) % :orient exact
                  constr_array{t,1} = [ 5, targt, obj_dir(targt) ];
               elseif (constr_type(1,t) == 6) % :orient range
                  arangex = orient_min_range + randi(orient_max_range-orient_min_range);
                  constr_array{t,1} = [ 6, targt, max(obj_dir(targt)-arangex,0), min(obj_dir(targt)+arangex,359) ];
               elseif (constr_type(1,t) == 7) % orient vector
                  arangex = orient_min_range + randi(orient_max_range-orient_min_range);
                  constr_array{t,1} = [ 6, targt, max(obj_dir(targt)-arangex,0), min(obj_dir(targt)+arangex,359) ];
               elseif (constr_type(1,t) == 8  ||  constr_type(1,t) == 10) % psi exact  or delta exact
                  if (issame == 1)   % delta exact
                     constr_array{t,1} = [ 10, targt, refern, ang1 ];
                  else
                     constr_array{t,1} = [ 8, targt, refern, ang1 ];
                  endif
               elseif (constr_type(1,t) == 9  ||  constr_type(1,t) == 11) % psi range  or delta range
                  arangex = orient_min_range + randi(orient_max_range-orient_min_range);
                  if (issame == 1)   % delta range
                     constr_array{t,1} = [ 11, targt, refern, max(ang1-arangex,0), min(ang1+arangex,359) ];
                  else
                     constr_array{t,1} = [ 9, targt, refern, max(ang1-arangex,0), min(ang1+arangex,359) ];
                  endif

               elseif (constr_type(1,t) == 12) % exact location
                  constr_array{t,1} = [ 12, targt, obj_locx(targt), obj_locy(targt) ];
               elseif (constr_type(1,t) == 13) % distance range
                  numdist = (obj_locx(targt)-obj_locx(refern))^2 + (obj_locy(targt)-obj_locy(refern))^2 ;
                  arangex = dist_range1 + randi(dist_range2-dist_range1);
                  constr_array{t,1} = [ 13, targt, refern, max(numdist-arangex,0), numdist+arangex ];
               elseif (constr_type(1,t) == 14) % exact distance
                  numdist = (obj_locx(targt)-obj_locx(refern))^2 + (obj_locy(targt)-obj_locy(refern))^2 ;
                  constr_array{t,1} = [ 14, targt, refern, numdist ];
               else
                  arangex = orient_min_range + randi(orient_max_range-orient_min_range);
                  constr_array{t,1} = [ 6, targt, max(obj_dir(targt)-arangex,0), min(obj_dir(targt)+arangex,359) ];
               endif

             endfor  % number of basic constraints

              res6 = writefile_asp(sample_basp_folder,nobj,nconstr,gridx,constr_array,0,0,zeros(nconstr,1),0,[],0,[]);
              pause(sleep_int);

                % Write explain constraints
              if (is_in_explain == 1) 
                 for explain_ind=1:length(explain_ratio)
                   nexplain = round(explain_ratio(explain_ind)*nconstr);
                   explain_inds = sort(randperm(nexplain,nexplain));   % random constraints to convert

                   trg_objexp = randi(nobj,nexplain); 
                   ref_objexp = randi(nobj,nexplain); 
                   econstr_type = randi(14,nexplain);  
                     % 1: same opra, 2: one-sided diff opra 3: diff opra 4: qualit dist, 5:orient exact, 6:orient range, 7:orient vector, 8.psi exact, 9:psi range, 10:delta exact, 11. delt range, 12:exact location, 13: distance range, 14: exact distance

                   explain_array = constr_array;   

                   for eind=1:nexplain
                      constr_ind = explain_inds(eind);
                      targt = trg_objexp(1,eind);
                      refern = ref_objexp(1,eind);
                      granul = 1 + randi(5);   %2*randi(3);

                      if (targt == refern) && (refern > 1)    % target and ref objects cannot be the same
                         refern = refern - 1;
                      elseif (targt == refern) && (refern == 1)
                         targt = 2;
                      endif

                      if (econstr_type(1,eind) == 1)   % same opra
                         explain_array{constr_ind,1} = [ 1, targt, refern, randi(4*granul)-1, granul ];
                      elseif (econstr_type(1,eind) == 2) % one sided diff opra
                         explain_array{constr_ind,1} = [ 2, targt, refern, randi(4*granul)-1, granul ];
                      elseif (econstr_type(1,eind) == 3) % diff opra
                         explain_array{constr_ind,1} = [ 3, targt, refern, randi(4*granul)-1, randi(4*granul)-1, granul ];
                      elseif (econstr_type(1,eind) == 4) % : qualit dist
                         explain_array{constr_ind,1} = [ 4, targt, refern, 1+randi(len_dist-1) ];
                      elseif (econstr_type(1,eind) == 5) % :orient exact
                         explain_array{constr_ind,1} = [ 5, targt, randi(360)-1 ];
                      elseif (econstr_type(1,eind) == 6) % :orient range
                         randangle = randi(360-orient_max_range,1,2) - 1;
                         arangex = orient_min_range + randi(orient_max_range-orient_min_range);
                         explain_array{constr_ind,1} = [ 6, targt, max(randangle-arangex,0), min(randangle+arangex,359) ];
                      elseif (econstr_type(1,eind) == 7) % orient vector
                         arangex = randi(nobj,1,2);
                         if (arangex(1) == arangex(2)) && (arangex(1) > 1)    % target and ref objects cannot be the same
                            arangex(1) = arangex(1) - 1;
                         elseif (arangex(1) == arangex(2)) && (arangex(1) == 1)
                            arangex(1) = 2;
                         endif
                         explain_array{constr_ind,1} = [ 7, targt, min(arangex), max(arangex) ];
                      elseif (econstr_type(1,eind) == 8) % psi exact
                         explain_array{constr_ind,1} = [ 8, targt, refern, randi(360)-1 ];
                      elseif (econstr_type(1,eind) == 9) % psi range
                         randangle = randi(360-orient_max_range,1,2) - 1;
                         arangex = orient_min_range + randi(orient_max_range-orient_min_range);
                         explain_array{constr_ind,1} = [ 9, targt, refern, max(randangle-arangex,0), min(randangle+arangex,359) ];
                      elseif (econstr_type(1,eind) == 10) % delta exact
                         explain_array{constr_ind,1} = [ 10, targt, refern, randi(360)-1 ];
                      elseif (econstr_type(1,eind) == 11) % delta range
                         randangle = randi(360-orient_max_range,1,2) - 1;
                         arangex = orient_min_range + randi(orient_max_range-orient_min_range);
                         explain_array{constr_ind,1} = [ 11, targt, refern, max(randangle-arangex,0), min(randangle+arangex,359) ];
                      elseif (econstr_type(1,eind) == 12) % exact location
                         explain_array{constr_ind,1} = [ 12, targt, randi(gridx), randi(gridx) ];
                      elseif (econstr_type(1,eind) == 13) % distance range
                         arangex = randi(2*gridx*gridx,1,2);
                         explain_array{constr_ind,1} = [ 13, targt, refern, min(arangex), max(arangex) ];
                      elseif (econstr_type(1,eind) == 14) % exact distance
                         explain_array{constr_ind,1} = [ 14, targt, refern, randi(2*gridx*gridx) ];
                      endif

                   endfor  % number of explain constraints

                   res71 = writefile_asp(sample_easp_folder,nobj,nconstr,gridx,explain_array,0,0,zeros(nconstr,1),0,[],nexplain,[]);
                   pause(sleep_int);

                 endfor
              endif

                % Write presumed (default) constraints
              if (is_in_def == 1) 
                 for def_ind=1:length(def_ratio)
                   ndef = round(def_ratio(def_ind)*nconstr);
                   def_inds = sort(randperm(nconstr,ndef));   % choose random constraints to convert to default

                   res41 = writefile_asp(sample_fasp_folder,nobj,nconstr,gridx,constr_array,0,0,zeros(nconstr,1),ndef,def_inds,0,[]);
                   pause(sleep_int);

                 endfor
              endif

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

                     disj_array = constr_array;   
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

                     res9 = writefile_asp(sample_dasp_folder,nobj,nconstr,gridx,disj_array,ndisj_constr,ndisjunc,constr_ind,0,[],0,orig_pair); %{rind,disj_ind}
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

  
  printf("------ Produced Consistent Instances----- \n");
  outcome =1;
  
  
endfunction
