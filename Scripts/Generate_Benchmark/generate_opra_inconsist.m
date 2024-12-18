function [outcome] = generate_opra_inconsist()



    function outcome = writefile_asp(afolder, nobject, nconstraint, carray, gr, num_disjconst, num_disjx, ind_array)

       if (num_disjconst == 0)
          fnamex = cstrcat(afolder,"/n", num2str(nobject), "_c", num2str(nconstraint), "_g", num2str(gr), "_basic_unsat.lp");
       else
          fnamex = cstrcat(afolder,"/n", num2str(nobject), "_c", num2str(nconstraint), "_g", num2str(gr), "_disj", num2str(num_disjconst), "x", num2str(num_disjx), "_unsat.lp");
       endif

       fptr = fopen(fnamex,"w");
       fprintf(fptr,"\n");
       fprintf(fptr,cstrcat("object(1..",num2str(nobject),").") );
       fprintf(fptr,"\n\n");
       fprintf(fptr,cstrcat("constraint_id(1..",num2str(nconstraint),").") );
       fprintf(fptr,"\n\n");

       %fprintf(fptr,cstrcat("num_digit(",num2str(dparam(1)),"). \n") );
       %fprintf(fptr,cstrcat("base(",num2str(dparam(2)),"). \n") );
       %fprintf(fptr,cstrcat("max_first_digit(",num2str(dparam(3)),"). \n\n") );

       % constr_array{t,1} = [ targt, refern, issame, sect1, sect2, grn, distind ];

       if (num_disjconst == 0)
          % all basic
          for cindx=1:nconstraint
             aconstr = carray{cindx,1};
             if (aconstr(3) == 1)
                   fprintf(fptr,cstrcat("same_opra(", num2str(cindx), ",", num2str(aconstr(2)), ",", num2str(aconstr(1)), ",", num2str(aconstr(4)), ",", num2str(aconstr(6)), "). \n" ) );
                else
                   fprintf(fptr,cstrcat("diff_opra(", num2str(cindx), ",", num2str(aconstr(1)), ",", num2str(aconstr(2)), ",", num2str(aconstr(4)), ",", num2str(aconstr(5)), ",", num2str(aconstr(6)), "). \n" ) );
                   fprintf(fptr,cstrcat("dist_rel(", num2str(cindx), ",", num2str(aconstr(1)), ",", num2str(aconstr(2)), ",", dist_name{aconstr(7)}, "). \n" ) );
             endif
             %fprintf(fptr,cstrcat("dist_rel(", num2str(aconstr(1)), ",", num2str(aconstr(2)), ",", dist_name{aconstr(7)}, "). \n" ) );
          endfor
       else          % not all basic
          dconstr_ind = 1;
          for cindx=1:nconstraint
             aconstr = carray{cindx,1};
             if(ind_array(cindx) == 0)       % basic
                if (aconstr(3) == 1)
                   fprintf(fptr,cstrcat("same_opra(", num2str(cindx), ",", num2str(aconstr(2)), ",", num2str(aconstr(1)), ",", num2str(aconstr(4)), ",", num2str(aconstr(6)), "). \n" ) );
                else
                   fprintf(fptr,cstrcat("diff_opra(", num2str(cindx), ",", num2str(aconstr(1)), ",", num2str(aconstr(2)), ",", num2str(aconstr(4)), ",", num2str(aconstr(5)), ",", num2str(aconstr(6)), "). \n" ) );
                   fprintf(fptr,cstrcat("dist_rel(", num2str(cindx), ",", num2str(aconstr(1)), ",", num2str(aconstr(2)), ",", dist_name{aconstr(7)}, "). \n" ) );
                endif
                %fprintf(fptr,cstrcat("dist_rel(", num2str(aconstr(1)), ",", num2str(aconstr(2)), ",", dist_name{aconstr(7)}, "). \n" ) );

             else              % disjunctive

                fprintf(fptr,cstrcat("disj_index(", num2str(cindx), ",1..", num2str(num_disjx), "). \n" ) );
                for (dindx=1:num_disjx)
                   %aconstr
                   if (aconstr{dindx}(3) == 0)       % diff opra
                      fprintf(fptr,cstrcat("disjrelation(", num2str(cindx), ",", num2str(dindx), ",", num2str(2), ",", num2str(aconstr{dindx}(1)), ",", num2str(aconstr{dindx}(2)), ",", num2str(aconstr{dindx}(4)), ",", num2str(aconstr{dindx}(5)), ",", num2str(aconstr{dindx}(6)), ",", dist_name{aconstr{dindx}(7)},  "). \n" ) );
                   elseif (aconstr{dindx}(3) == 1)   % same opra
                      fprintf(fptr,cstrcat("disjrelation(", num2str(cindx), ",", num2str(dindx), ",", num2str(1), ",", num2str(aconstr{dindx}(1)), ",", num2str(aconstr{dindx}(2)), ",", num2str(aconstr{dindx}(4)), ",", num2str(aconstr{dindx}(6)), "). \n" ) ); 
                   endif
                endfor
                dconstr_ind = dconstr_ind+1;
             endif
          endfor
       endif

       fclose(fptr);
       outcome = 1;
       
    endfunction





  dist_name = {"on", "verynear", "near", "medium", "far", "veryfar"};
  
  startindex = 1;
  num_samples = 100;

  % Number of objects, constraints and Granul
  num_objects = [10, 20, 40, 60];  %, 50, 100, 200];  %, 50, 100, 200, 400, 500];
  num_constr = {[10, 20], [20, 40], [40, 60], [60, 80]};
  num_granul = { {[2,4,6],[2,4]}, {[2,4,6],[2,4]}, {[2,4],[2,4]}, {[2,4],[2,4]} };

  disj_objects = [10, 20, 40];   %, 100, 200];  %, 50, 100, 200, 400, 500];
  disj_constr = {[10, 20], [20], [40]};
  disj_granul = {[2, 4], [2, 4], [2]}; 

  num_disjunc = [2,4,8];
  disjunc_ratio =  [0.2, 0.4];  % [0.1, 0.2];
  max_disjunc = 8;


  % Distance relations 
  distrel = 1:6;    % 0:5
  
  len_dist = length(distrel);

  
  aspfolder = "./ASP";

  consist_folder = "/Consistent";
  inconsist_folder = "/Inconsistent";

  sample_folder = "/Sample";

  basic_folder = "/Basic";
  disj_folder = "/Disjunctive";

  sleep_int = 2.0;
  sleep_sample = 11.0;
  
  len_objects = length(num_objects);
  len_disj_objects = length(disj_objects);
  len_disjunc = length(num_disjunc);
  len_ratio = length(disjunc_ratio);

   

  printf("------ RANDOM SAMPLES----- \n");
  printf("Finished Sample ");
  
  for m=1:num_samples
     sample_asp_folder = cstrcat(aspfolder,inconsist_folder,sample_folder,num2str(m+startindex-1));

     res1a = mkdir(sample_asp_folder);
     pause(sleep_int);

     sample_basp_folder = cstrcat(sample_asp_folder,basic_folder);
     sample_dasp_folder = cstrcat(sample_asp_folder,disj_folder);

     res3a = mkdir(sample_basp_folder);
     res3b = mkdir(sample_dasp_folder);
     pause(sleep_int);
     
     for n=1:len_objects
        nobj = num_objects(n); 
        len_constr = length(num_constr{n});

        for k=1:len_constr
           nconstr = num_constr{n}(k);  
           constr_granul = (num_granul{n}){k};
           len_granul = length(constr_granul);       
      
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

           if (is_in_disj == 1)  
              constr_ind = cell(len_ratio,len_disjunc);
              disj_constr_ind = cell(len_ratio,len_disjunc);
              %target_objs = cell(len_ratio,len_disjunc);
              %refer_objs = cell(len_ratio,len_disjunc);
              rand_same = cell(len_ratio,len_disjunc);
              rangle1 = cell(len_ratio,len_disjunc);
              rangle2 = cell(len_ratio,len_disjunc);
              rand_dist = cell(len_ratio,len_disjunc);
 
             for rind=1:len_ratio
                ndisj_constr = round(disjunc_ratio(rind)*nconstr);

                for disj_ind=1:len_disjunc
                   ndisjunc = num_disjunc(disj_ind);

                   constr_ind{rind,disj_ind} = zeros(nconstr,1);   % 0 for basic
                   disj_constr_ind{rind,disj_ind} = sort(randperm(nconstr,ndisj_constr));   % choose random constraints to convert
                   (constr_ind{rind,disj_ind})(disj_constr_ind{rind,disj_ind}) = 1;   % choose random constraints to convert to disjunctive

                   rand_same{rind,disj_ind} = randi(2,ndisj_constr,ndisjunc-1)-1;
                   rangle1{rind,disj_ind} = randi(360,ndisj_constr,ndisjunc-1) - 1;
                   rangle2{rind,disj_ind} = randi(360,ndisj_constr,ndisjunc-1) - 1;
                   rand_dist{rind,disj_ind} = 1+randi(len_dist-1,ndisj_constr,ndisjunc-1);

                 endfor    % for each disjunct
              endfor    % for each disjunctive constraint
           endif

           % create objects
              % generate random angles and is_same
           rand_s = randi(2,1,nconstr) - 1;
           rangle3 = randi(360,1,nconstr) - 1;
           rangle4 = randi(360,1,nconstr) - 1;
           rdist = 1+randi(len_dist-1,1,nconstr);
           obj_dir = randi(360,1,nobj) - 1; 

           apermut = randperm(nobj*(nobj-1), nconstr);
           for tind=1:nconstr
                trg_obj(tind) = 1 + floor((apermut(tind)-1) / (nobj-1));  % x index row id ranges from 1 to nobj, 1 is topmost
                ref_obj(tind) = 1 + mod(apermut(tind)-1, nobj-1);    % yindex column id 1 is leftmost
                if (ref_obj(tind) >= trg_obj(tind))  
                   ref_obj(tind) = ref_obj(tind) + 1;   % escape diagonal
                endif
           endfor


           for z=1:len_granul
              granul = constr_granul(z);
              quadr = 180 / granul;
              %granul = ((num_granul{n}){k})(z);      

              constr_array = cell(nconstr,1); 
              for t=1:nconstr
                 targt = trg_obj(t);
                 refern = ref_obj(t);

                 sect_coef1 = floor(rangle3(t)/quadr);
                 sect_mod1 = mod(rangle3(t), quadr);
      
                 if (sect_mod1 > 0) 
                    sect1 = 2*sect_coef1 + 1;
                 else
                    sect1 = 2*sect_coef1;
                 endif

                 sect_coef2 = floor(rangle4(t)/quadr);
                 sect_mod2 = mod(rangle4(t), quadr);
         
                 if (sect_mod2 > 0) 
                    sect2 = 2*sect_coef2 + 1;
                 else
                    sect2 = 2*sect_coef2;
                 endif

                 constr_array{t,1} = [ targt, refern, rand_s(t), sect1, sect2, granul, rdist(t) ];
  
              endfor  % number of basic constraints

              res6 = writefile_asp(sample_basp_folder,nobj,nconstr,constr_array,granul,0,0,zeros(nconstr,1));
              pause(sleep_int);


              if (is_in_disj == 1) 
                if ismember(granul,disj_granul{nindex})
 
                 for rind=1:len_ratio
                   ndisj_constr = round(disjunc_ratio(rind)*nconstr);

                   for disj_ind=1:len_disjunc
                      ndisjunc = num_disjunc(disj_ind);

                      disj_array = constr_array;   

                         % create disjunctive constraints
                      for cind=1:ndisj_constr
                           % choose constraints to convert to disjunctive
                         constr_id = (disj_constr_ind{rind,disj_ind})(cind);
                         disj_array{constr_id,1} = cell(ndisjunc,1);
                         disj_array{constr_id,1}{1} = constr_array{constr_id,1};
                         for dindex=1:(ndisjunc-1)
                             trgd = trg_obj(constr_id);    
                             refd = ref_obj(constr_id);     

                             sect_coef1 = floor((rangle1{rind,disj_ind})(cind,dindex)/quadr);
                             sect_mod1 = mod((rangle1{rind,disj_ind})(cind,dindex), quadr);
                             if (sect_mod1 > 0) 
                                sect1 = 2*sect_coef1 + 1;
                             else
                                sect1 = 2*sect_coef1;
                             endif

                             sect_coef2 = floor((rangle2{rind,disj_ind})(cind,dindex)/quadr);
                             sect_mod2 = mod((rangle2{rind,disj_ind})(cind,dindex), quadr);
                             if (sect_mod2 > 0) 
                                sect2 = 2*sect_coef2 + 1;
                             else
                                sect2 = 2*sect_coef2;
                             endif
 
                             disj_array{constr_id,1}{dindex+1} = [trgd, refd, (rand_same{rind,disj_ind})(cind,dindex), sect1, sect2, granul, (rand_dist{rind,disj_ind})(cind,dindex)]; 
                         endfor    % for each disjunct
                      endfor    % for each disjunctive constraint

                      res9 = writefile_asp(sample_dasp_folder,nobj,nconstr,disj_array,granul,ndisj_constr,ndisjunc,constr_ind{rind,disj_ind});
                      pause(sleep_int);

                   endfor    % num disjuncs
                 endfor   %  disjunctive ratio

              endif  % if in disj granul
             endif  % if in disjunc size

          endfor    % num granul

       endfor   %  num constraints

     endfor       % num objects

     printf( cstrcat(num2str(m+startindex-1), " ") );            %printf(cstrcat("Finished Sample ", num2str(m+startindex-1)) );
     pause(sleep_sample);

   endfor       % samples 

  
  printf("------ Produced Inconsistent Instances----- \n");
  outcome =1;
  
endfunction
