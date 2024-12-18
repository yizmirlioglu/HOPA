
import os
import subprocess
import time
import threading
import psutil
import itertools
import random
import string
import copy
import math
import statistics


DEBUG = True


    
def test_inst(solvername, instfile, cspfile, solvertype):

   global continue_measure
   global used_mem_array
   global init_measure_phase
   global is_intermiss
   global cur_process

   time.sleep(interval)
   init_mem = used_mem_array[-1]
   used_mem_array = used_mem_array[-last_ind:-1] + [used_mem_array[-1]]
   init_measure_phase = True

   start = time.time()

   if (solvertype == 1):   # flatzinc
      comp_command = compile_command + " " + instfile + ".mzn " + cspfile  # + ".mzn"
      p1 = subprocess.Popen(comp_command, stdout=subprocess.PIPE, shell=True)
      (outputx1, err1) = p1.communicate()
      p1.wait()
      outputx1 = str(outputx1)

      command = solvername + " " + instfile + ".fzn" # + " -s --time-limit " + str(timeout_value)
   elif (solvertype == 0):   
      command = solvers_upper_folder + "/" + solvername + " -s --time-limit " + str(timeout_value*1000) + " " + cspfile + " " + instfile + ".mzn " 

   p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

   (outputx, err) = p.communicate()
   p.wait()

   outputx = str(outputx)

   ctime = time.time() - start
   ctimefx = float(format(ctime, '.2f'))

   final_mem = max(used_mem_array)
   mem_consume = final_mem - init_mem

   cpu_time = -1

   unsatindex = outputx.find("UNSATISFIABLE")
   unsatindex2 = outputx.find("Unsatisfiable")
   modelindex = outputx.find("model inconsistency detected")

   is_model = 0
   if (modelindex >= 0):
      is_model = 1

   is_tout = 0
   if (ctime >= (timeout_value-2)):
      is_tout = 1

   is_sat = 1
   if (unsatindex >= 0) or (unsatindex2 >= 0) or (modelindex >= 0) or (is_tout == 1):
      is_sat = 0


   return ctimefx, cpu_time, is_sat, is_model, is_tout, mem_consume





def extract_values(aspout):

   tsind = 0   #toutput.rfind("SATISFIABLE: ")
   timeindx = aspout.find("Time   ", tsind)
   sp_ind1 = aspout.find(":", timeindx+4)
   sp_ind2 = aspout.find("s", sp_ind1+1)
   total_text_all = aspout[(sp_ind1+1):sp_ind2]
   total_float = float(total_text_all)
   total_time = float(format(total_float, '.3f'))

   search_ind = aspout.find("Solving:", timeindx)
   sp_ind = aspout.find("s", search_ind+8)
   search_text_all = aspout[(search_ind+9):sp_ind]
   search_float = float(search_text_all)
   search_time = float(format(search_float, '.3f'))

   cpu_ind = aspout.find("CPU Time", timeindx)
   sp_ind1 = aspout.find(":", cpu_ind+8)
   sp_ind2 = aspout.find("s", sp_ind1+1)
   cpu_text_all = aspout[(sp_ind1+1):sp_ind2]
   cpu_float = float(cpu_text_all)
   cpu_time = float(format(cpu_float, '.3f'))
   #sp_ind = aspout.find("s", cpu_ind+8)
   #cpu_text_all = aspout[(cpu_ind+9):sp_ind]

   atoms_ind = aspout.find("Atoms")
   atoms_text = aspout[(atoms_ind+14):(atoms_ind+23)]
   atoms_count = int(atoms_text)

   rules_ind = aspout.find("Rules")
   rules_text = aspout[(rules_ind+14):(rules_ind+23)]
   rules_count = int(rules_text)
   vars_ind = aspout.find("Variables")
   vars_text = aspout[(vars_ind+14):(vars_ind+23)]
   vars_count = int(vars_text)
   constr_ind = aspout.find("Constraints")
   constr_text = aspout[(constr_ind+14):(constr_ind+23)]
   constr_count = int(constr_text)

   ground_time = total_time - search_time

   return total_time, cpu_time, search_time, ground_time, atoms_count, rules_count, vars_count, constr_count




def measure_memory():

   global continue_measure
   global is_intermiss
   global used_mem_array
   global init_measure_phase
   global cur_process #= psutil.Process(os.getpid())

   phase_counter = 1
   #used_mem_array = []

   while(continue_measure):

      if (not is_intermiss):
         #used_mb=int(cur_process.memory_info().rss / 1048576.0)   #  1024.0*1024.0
         mem_info = psutil.virtual_memory()
         used_mb = float(format(mem_info.used / 1073741824.0, '.6f'))   # mb= 1048576.0)   #  1024.0*1024.0

         used_mem_array.append(used_mb)

         if (init_measure_phase and phase_counter < len_init_phase):
             time.sleep( init_wait_periods[phase_counter-1] )
             phase_counter = phase_counter + 1
             #print("Init mem phase")
         elif (init_measure_phase and phase_counter == len_init_phase):
             init_measure_phase = False
             time.sleep( init_wait_periods[phase_counter-1] )
             phase_counter = 1
             #print("Exiting init mem phase")
         else:
             #print(used_mem_array)
             #print("Regular mem phase")
             time.sleep(mem_step)    #intermiss_mem_step)
      else:
         print("Memory measure entered intermission")
         time.sleep( intermiss+intermiss_mem_step )    # mem_step,   intermiss+0.6
         init_measure_phase = False
         phase_counter = 1




print("Real Test has started ... ")


# Parameters for measuring memory (RAM)
mem_step = 0.5            # sample period for memory
intermiss_mem_step = 1.0  # wait period during intermission for memory measure
refr = 50           # empty mem array after this length
last_ind = 2        # leave last elements after emptying


global used_mem_array
global init_measure_phase
global continue_measure
global is_intermiss
global cur_process

cur_process = psutil.Process(os.getpid())

used_mem_array = []   # memory used (global var)
continue_measure = True
init_measure_phase = False
is_intermiss = False

init_wait_periods = [0.05, 0.05, 0.1]
len_init_phase = len(init_wait_periods)

mem_infox = psutil.virtual_memory()
used_mem_mb = float( format( mem_infox.used / 1073741824.0, '.6f' ) )   # mb= 1048576.0)   #  1024.0*1024.0  

used_mem_array.append(used_mem_mb)


# start memory thread

t1 = threading.Thread(target=measure_memory, args=())   #  args=(10,)
t1.start()



# INCONSISTENT TEST STARTS HERE

sleep_int = 10.0
sleep_sample = 60.0
timeout_value = 100
intermiss = 7000.0

startsample = 1
nsample = 100

num_objects = [10, 20, 40, 60];   #
num_constr = [ [10,20], [20,40], [40], [60,80] ]  
all_granuls = [ [[4,6],[4]], [[2,4],[2]], [[2]], [[2],[2]] ]

all_grids =  [ ["grid100.lp ", "grid50.lp "], ["grid50.lp "], ["grid50.lp "], ["grid25.lp "] ]   

all_resol =  [ ["deg2.lp ", "deg3.lp "], ["deg3.lp "], ["deg3.lp "], ["deg3.lp "] ]   

grid_names =  [ ["g100 deg2", "g50 deg3"], ["g50 deg3"], ["g50 deg3"], ["g25 deg3"] ]  



engines = ["./clingo "]
solver_files = ["opra_consistency.lp"]    # "main_canonical.lp disjunctive1.lp guide_cdc.lp", "main_full.lp disjunctive1.lp guide_cdc.lp"


aspfolder = "./ASP"
cspfolder = "./MiniZinc"

consist_folder = "Consistent"
inconsist_folder = "Inconsistent"

sample_folder = "Sample"

basic_folder = "Basic"
disj_folder = "Disjunctive"
infer_folder = "Infer"

continue_measure = True
init_measure_phase = False
is_intermiss = False


instances_upper_folder = "./OPRA_Instances"  #"./Consistent"

solver_files_folder = "./ASP"


result_folder = os.path.join("Results")

rfname_unsat = "Basic_OPRA_Inconsist_Results.csv"
result_file_unsat = os.path.join(result_folder, rfname_unsat)


time_sample_results = []
sat_sample_results = []
timeout_sample_results = []

search_time_results = []
cpu_time_results = []

memory_results = []
atoms_results = []
rules_results = []
vars_results = []
constr_results = []


engineind = 0
for aengine in engines:

   time_sample_results.append([])
   sat_sample_results.append([])
   timeout_sample_results.append([])
   search_time_results.append([])
   cpu_time_results.append([])
   memory_results.append([])
   atoms_results.append([])
   rules_results.append([])
   vars_results.append([])
   constr_results.append([])

   solverind = 0
   for asolver in solver_files:
      #stype = solver_type[solverind]

      (time_sample_results[engineind]).append([])
      (sat_sample_results[engineind]).append([])
      #(gridx_sample_results[engineind]).append([])
      #(gridy_sample_results[engineind]).append([])
      (timeout_sample_results[engineind]).append([])
      #(grid_time_results[engineind]).append([])
      (search_time_results[engineind]).append([])
      (cpu_time_results[engineind]).append([])
      (memory_results[engineind]).append([])
      (atoms_results[engineind]).append([])
      (rules_results[engineind]).append([])
      (vars_results[engineind]).append([])
      (constr_results[engineind]).append([])

      objind = 0
      for nobj in num_objects:
         constr_list = num_constr[objind]
         grids = all_grids[objind]
         resols = all_resol[objind]

         ((time_sample_results[engineind])[solverind]).append([])
         ((sat_sample_results[engineind])[solverind]).append([])
         #((gridx_sample_results[engineind])[solverind]).append([])
         #((gridy_sample_results[engineind])[solverind]).append([])
         ((timeout_sample_results[engineind])[solverind]).append([])
         #((grid_time_results[engineind])[solverind]).append([])
         ((search_time_results[engineind])[solverind]).append([])
         ((cpu_time_results[engineind])[solverind]).append([])
         ((memory_results[engineind])[solverind]).append([])
         ((atoms_results[engineind])[solverind]).append([])
         ((rules_results[engineind])[solverind]).append([])
         ((vars_results[engineind])[solverind]).append([])
         ((constr_results[engineind])[solverind]).append([])

         cind = 0
         for aconstr in constr_list:
            granuls = (all_granuls[objind])[cind]

            (((time_sample_results[engineind])[solverind])[objind]).append([])
            (((sat_sample_results[engineind])[solverind])[objind]).append([])
            #(((gridx_sample_results[engineind])[solverind])[objind]).append([])
            #(((gridy_sample_results[engineind])[solverind])[objind]).append([])
            (((timeout_sample_results[engineind])[solverind])[objind]).append([])
            #(((grid_time_results[engineind])[solverind])[objind]).append([])
            (((search_time_results[engineind])[solverind])[objind]).append([])
            (((cpu_time_results[engineind])[solverind])[objind]).append([])
            (((memory_results[engineind])[solverind])[objind]).append([])
            (((atoms_results[engineind])[solverind])[objind]).append([])
            (((rules_results[engineind])[solverind])[objind]).append([])
            (((vars_results[engineind])[solverind])[objind]).append([])
            (((constr_results[engineind])[solverind])[objind]).append([])

            gind = 0
            for agranul in granuls:
               ((((time_sample_results[engineind])[solverind])[objind])[cind]).append([])
               ((((sat_sample_results[engineind])[solverind])[objind])[cind]).append([])
               #((((gridx_sample_results[engineind])[solverind])[objind])[cind]).append([])
               #((((gridy_sample_results[engineind])[solverind])[objind])[cind]).append([])
               ((((timeout_sample_results[engineind])[solverind])[objind])[cind]).append([])
               #((((grid_time_results[engineind])[solverind])[objind])[cind]).append([])
               ((((search_time_results[engineind])[solverind])[objind])[cind]).append([])
               ((((cpu_time_results[engineind])[solverind])[objind])[cind]).append([])
               ((((memory_results[engineind])[solverind])[objind])[cind]).append([])
               ((((atoms_results[engineind])[solverind])[objind])[cind]).append([])
               ((((rules_results[engineind])[solverind])[objind])[cind]).append([])
               ((((vars_results[engineind])[solverind])[objind])[cind]).append([])
               ((((constr_results[engineind])[solverind])[objind])[cind]).append([])

               #num_disjconstr = int(agranul*aconstr)

               gridind = 0
               for agrid in grids:
                  (((((time_sample_results[engineind])[solverind])[objind])[cind])[gind]).append([])
                  (((((sat_sample_results[engineind])[solverind])[objind])[cind])[gind]).append([])
                  #(((((gridx_sample_results[engineind])[solverind])[objind])[cind])[gind]).append([])
                  #(((((gridy_sample_results[engineind])[solverind])[objind])[cind])[gind]).append([])
                  (((((timeout_sample_results[engineind])[solverind])[objind])[cind])[gind]).append([])
                  #(((((grid_time_results[engineind])[solverind])[objind])[cind])[gind]).append([])
                  (((((search_time_results[engineind])[solverind])[objind])[cind])[gind]).append([])
                  (((((cpu_time_results[engineind])[solverind])[objind])[cind])[gind]).append([])
                  (((((memory_results[engineind])[solverind])[objind])[cind])[gind]).append([])
                  (((((atoms_results[engineind])[solverind])[objind])[cind])[gind]).append([])
                  (((((rules_results[engineind])[solverind])[objind])[cind])[gind]).append([])
                  (((((vars_results[engineind])[solverind])[objind])[cind])[gind]).append([])
                  (((((constr_results[engineind])[solverind])[objind])[cind])[gind]).append([])

                  aresol = resols[gridind]
                  inst_name_unsat = "n" + str(nobj) + "_c" + str(aconstr) + "_g" + str(agranul) + "_basic_unsat.lp"

                  sampleind = 0
                  while sampleind < nsample:

                     instance_folder = os.path.join(instances_upper_folder, inconsist_folder, sample_folder+str(sampleind+startsample), basic_folder )

                     asolver_file = solver_files_folder + "/" + asolver    #os.path.join(solver_files_folder, asolver)

                     file_inst = os.path.join(instance_folder, inst_name_unsat)

                     tcommand = aengine + " " + asolver_file + " " +  solver_files_folder + "/" + str(agrid) + solver_files_folder + "/" + str(aresol) + solver_files_folder + "/" +  "tangent.lp " + file_inst + " --time-limit=" + str(timeout_value) + " --warn none --stats=2"


                     begintime = time.time()

                     # Test instance now
                     init_mem = used_mem_array[-1]
                     used_mem_array = used_mem_array[-last_ind:-1] + [used_mem_array[-1]]
                     init_measure_phase = True

                     p = subprocess.Popen(tcommand, stdout=subprocess.PIPE, shell=True)
                     (output, err) = p.communicate()
                     p.wait()
                     ctime = time.time() - begintime
                     ctimefx = float(format(ctime, '.3f'))
                     #solvetime = float(format(ctime - gridtime, '.3f'))
                     #gridfx = float(format(gridtime, '.3f'))

                     final_mem = max(used_mem_array)
                     mem_consume = max(final_mem - init_mem, 0)

                     toutput = str(output)
                     #print(toutput)

                     totalt, cput, searcht, grndt, atomc, rulec, varc, constrc = extract_values(toutput)

                     unsatindex = toutput.find("UNSATISFIABLE")

                     if (ctime < (timeout_value-1.2)):  #  no timeout
                        ((((((time_sample_results[engineind])[solverind])[objind])[cind])[gind])[gridind]).append(totalt)
                        #((((((grid_time_results[engineind])[solverind])[objind])[cind])[gind])[gridind]).append(gridfx)
                        ((((((search_time_results[engineind])[solverind])[objind])[cind])[gind])[gridind]).append(searcht)
                        ((((((cpu_time_results[engineind])[solverind])[objind])[cind])[gind])[gridind]).append(cput)  # solvetime
                        ((((((timeout_sample_results[engineind])[solverind])[objind])[cind])[gind])[gridind]).append(0) 

                        if (unsatindex >= 0):
                           ((((((sat_sample_results[engineind])[solverind])[objind])[cind])[gind])[gridind]).append(0)
                        else:
                           ((((((sat_sample_results[engineind])[solverind])[objind])[cind])[gind])[gridind]).append(1)

                     else:
                        ((((((timeout_sample_results[engineind])[solverind])[objind])[cind])[gind])[gridind]).append(1) 


                     ((((((memory_results[engineind])[solverind])[objind])[cind])[gind])[gridind]).append(mem_consume)
                     ((((((atoms_results[engineind])[solverind])[objind])[cind])[gind])[gridind]).append(atomc)
                     ((((((rules_results[engineind])[solverind])[objind])[cind])[gind])[gridind]).append(rulec)
                     ((((((vars_results[engineind])[solverind])[objind])[cind])[gind])[gridind]).append(varc)
                     ((((((constr_results[engineind])[solverind])[objind])[cind])[gind])[gridind]).append(constrc)
                     #((((((gridx_sample_results[engineind])[solverind])[objind])[cind])[gind])[gridind]).append(gridx)
                     #((((((gridy_sample_results[engineind])[solverind])[objind])[cind])[gind])[gridind]).append(gridy)


                     time.sleep(sleep_int)
                     sampleind = sampleind + 1

                  time.sleep(sleep_sample)
                  gridind = gridind + 1
               gind = gind + 1

            print("Finished all samples for Object and Constraint: ", str(nobj), " ", str(aconstr), " ")

            cind = cind + 1

         print("Finished Num object: ", str(nobj))
         objind = objind + 1

      print("Finished Solver File ID", str(solverind+1))

      solverind = solverind + 1

   engineind = engineind + 1


#  WRITE RESULTS


print("WRITING RESULTS ------", "\n")

columnslot = 32    # 6*len(solvers) + 6

excelfile1 = open(result_file_unsat, "w")

slotind = 0
while slotind < (columnslot-1):
    excelfile1.write(",")
    slotind = slotind + 1  

excelfile1.write("\n")

slotind = 0
while slotind < (columnslot-1):
    excelfile1.write(",")
    slotind = slotind + 1 
excelfile1.write("\n")

engineind = 0
for aengine in engines:

   solverind = 0
   for asolver in solver_files:

      excelfile1.write(",,,," + str(aengine) + ",,,,,,,,,,,,,,,,,,,,,")
      excelfile1.write("\n")
      excelfile1.write(",,,," + str(asolver) + ",,,,,,,,,,,,,,,,,,,,,")
      excelfile1.write("\n\n")

      excelfile1.write("N,Constr,Granul,Grid,,Total,Ground,Search,Median,StdDev,MaxTime,MinTime,Satisfiable,Timeout,,Peak Memory,Aver Memory,Median Memory,,Peak Atoms,Peak Rules,Peak Vars,Peak Constr,,Median Atoms,Median Rules,Median Vars,Median Constr,,")
      excelfile1.write("\n")

      objind = 0
      for nobj in num_objects:
         constr_list = num_constr[objind]
         #granuls = all_granuls[objind]
         grids = all_grids[objind]
         grid_cases = grid_names[objind]

         cind = 0
         for aconstr in constr_list:
            granuls = (all_granuls[objind])[cind]

            gind = 0
            for agranul in granuls:

               gridind = 0
               for agrid in grid_cases:

                  if len( (((((time_sample_results[engineind])[solverind])[objind])[cind])[gind])[gridind] ) > 1:
                     avertime = float(format( statistics.mean( (((((time_sample_results[engineind])[solverind])[objind])[cind])[gind])[gridind] ), '.3f') )
                     medx = float(format( statistics.median( (((((time_sample_results[engineind])[solverind])[objind])[cind])[gind])[gridind] ), '.3f') )
                     stdx = float(format( statistics.stdev( (((((time_sample_results[engineind])[solverind])[objind])[cind])[gind])[gridind] ), '.3f') )
                     maxtime = float(format( max( (((((time_sample_results[engineind])[solverind])[objind])[cind])[gind])[gridind] ), '.3f') )
                     mintime = float(format( min( (((((time_sample_results[engineind])[solverind])[objind])[cind])[gind])[gridind] ), '.3f') )

                     avercpu = float(format( statistics.mean( (((((cpu_time_results[engineind])[solverind])[objind])[cind])[gind])[gridind] ), '.3f') )
                     aversearch = float(format( statistics.mean( (((((search_time_results[engineind])[solverind])[objind])[cind])[gind])[gridind] ), '.3f') )
                     averground = avertime - aversearch
                  else:
                     avertime = -1
                     medx = -1
                     stdx = -1
                     maxtime = -1
                     mintime = -1
                     avercpu = -1
                     aversearch = -1
                     averground = -1

                  satsum = sum( (((((sat_sample_results[engineind])[solverind])[objind])[cind])[gind])[gridind] )
                  #gridx_val = float(format( statistics.median( ((((((gridx_sample_results[engineind])[solverind])[objind])[cind])[gind])[gridind]) ), '.1f') )
                  #gridy_val = float(format( statistics.median( ((((((gridy_sample_results[engineind])[solverind])[objind])[cind])[gind])[gridind]) ), '.1f') )
                  timeoutsum = sum( (((((timeout_sample_results[engineind])[solverind])[objind])[cind])[gind])[gridind] )

                  peak_mem = float(format( max( (((((memory_results[engineind])[solverind])[objind])[cind])[gind])[gridind] ), '.3f') )
                  aver_mem = float(format( statistics.mean( (((((memory_results[engineind])[solverind])[objind])[cind])[gind])[gridind] ), '.3f') )
                  median_mem = float(format( statistics.median( (((((memory_results[engineind])[solverind])[objind])[cind])[gind])[gridind] ), '.3f') )

                  peak_atoms = max( (((((atoms_results[engineind])[solverind])[objind])[cind])[gind])[gridind] )
                  peak_rules = max( (((((rules_results[engineind])[solverind])[objind])[cind])[gind])[gridind] )
                  peak_vars = max( (((((vars_results[engineind])[solverind])[objind])[cind])[gind])[gridind] )
                  peak_constr = max( (((((constr_results[engineind])[solverind])[objind])[cind])[gind])[gridind] )

                  median_atoms = statistics.median_high( (((((atoms_results[engineind])[solverind])[objind])[cind])[gind])[gridind] )
                  median_rules = statistics.median_high( (((((rules_results[engineind])[solverind])[objind])[cind])[gind])[gridind] )
                  median_vars = statistics.median_high( (((((vars_results[engineind])[solverind])[objind])[cind])[gind])[gridind] )
                  median_constr = statistics.median_high( (((((constr_results[engineind])[solverind])[objind])[cind])[gind])[gridind] )

                  excelfile1.write(str(nobj) + "," + str(aconstr) + "," + str(agranul) + "," + str(agrid) + ",," + str(avertime) + "," + str(averground) + "," + str(aversearch) + "," + str(medx) + "," + str(stdx) + ","  + str(maxtime) + "," + str(mintime) + "," + str(satsum) + "," + str(timeoutsum) + ",," + str(peak_mem) + "," + str(aver_mem) + "," + str(median_mem) + ",," + str(peak_atoms) + "," + str(peak_rules) + "," + str(peak_vars) + "," + str(peak_constr) + ",," + str(median_atoms) + "," + str(median_rules) + "," + str(median_vars) + "," + str(median_constr) + ",," )

                  excelfile1.write("\n")

                  gridind = gridind + 1
               gind = gind + 1

            cind = cind +1
         objind = objind + 1

      excelfile1.write("\n")
      #time.sleep(sleep_int)

      slotind = 0
      while slotind < (columnslot-1):
          excelfile1.write(",")
          slotind = slotind + 1 

      excelfile1.write("\n")

      slotind = 0
      while slotind < (columnslot-1):
          excelfile1.write(",")
          slotind = slotind + 1 
      excelfile1.write("\n")

      solverind = solverind + 1

   engineind = engineind + 1

excelfile1.close()





continue_measure = False         # stop measuring memory
time.sleep(sleep_int+mem_step)

t1.join()       # wait until thread 1 is completely executed

time.sleep(sleep_int)

