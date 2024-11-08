[1mdiff --git a/projection_method/parsers/gtensor_calculation.py b/projection_method/parsers/gtensor_calculation.py[m
[1mindex f1aa307..bcdf2b7 100644[m
[1m--- a/projection_method/parsers/gtensor_calculation.py[m
[1m+++ b/projection_method/parsers/gtensor_calculation.py[m
[36m@@ -19,20 +19,20 @@[m [mfile = str(sys.argv[1])[m
 [m
 ######## G-TENSOR CALCULATION ########[m
 g_calculation = 1[m
[31m-ppm = 1 # 0: ppt; 1: ppm[m
[32m+[m[32mppm = 0 # 0: ppt; 1: ppm[m
 [m
 ######## G-TENSOR ANALYSIS ########[m
 excitanalysis_gvalue_cut = 0 # =0: not calculate; ≠0: cut-off between ground-excited states (% of maximum g-value in each dim)[m
 [m
 ######## EXCITED STATES ANALYSIS ########[m
 excitanalysis = 0[m
[31m-cutoffamp = 0.9 # cut-off for configurations amplitude[m
[31m-excitanalysis_plot = 1[m
[32m+[m[32mcutoffamp = 0 # cut-off for configurations amplitude[m
[32m+[m[32mexcitanalysis_plot = 0[m
 [m
 ######## SOS PLOTS ########[m
 sos_analysis = 0 # SOS g-tensor plot: g-tensor calculation with n states[m
[31m-amp_cutoff = 0.1[m
[31m-g_estimation = 1[m
[32m+[m[32mamp_cutoff = 0[m
[32m+[m[32mg_estimation = 0[m
 [m
 [m
 def extract_data_from_json(filee):[m
[36m@@ -347,19 +347,19 @@[m [mdef print_g_calculation(filee, totalstates, gtensor, ppms, spin_list):[m
                 multip_dict[new_key] = counter_dict[old_key][m
         return multip_dict[m
     [m
[31m-    multipdict = count_spin_states(spin_list)[m
[31m-    print("---------------------")[m
[31m-    print(" G-SHIFT RESULTS")[m
[31m-    print("---------------------")[m
[31m-    print("File selected: ", filee)[m
[31m-    print("Number of states: ", totalstates)[m
[31m-    print("States multiplicity: ")[m
[31m-    [print(" ", key, value, end=',') for key, value in multipdict.items()][m
[31m-    print()[m
[31m-    if ppms == 0:[m
[31m-        print('g-factor (x y z dimensions) in ppt:')[m
[31m-    elif ppms == 1:[m
[31m-        print('g-factor (x y z dimensions) in ppm:')[m
[32m+[m[32m    # multipdict = count_spin_states(spin_list)[m
[32m+[m[32m    # print("---------------------")[m
[32m+[m[32m    # print(" G-SHIFT RESULTS")[m
[32m+[m[32m    # print("---------------------")[m
[32m+[m[32m    # print("File selected: ", filee)[m
[32m+[m[32m    # print("Number of states: ", totalstates)[m
[32m+[m[32m    # print("States multiplicity: ")[m
[32m+[m[32m    # [print(" ", key, value, end=',') for key, value in multipdict.items()][m
[32m+[m[32m    # print()[m
[32m+[m[32m    # if ppms == 0:[m
[32m+[m[32m    #     print('g-factor (x y z dimensions) in ppt:')[m
[32m+[m[32m    # elif ppms == 1:[m
[32m+[m[32m    #     print('g-factor (x y z dimensions) in ppm:')[m
     print(np.round(gtensor[0].real, 3), np.round(gtensor[1].real, 3),[m
           np.round(gtensor[2].real, 3))[m
     print('')[m
[36m@@ -673,10 +673,10 @@[m [mif g_calculation == 1:[m
 [m
     print_g_calculation(file, nstates, gshift, ppm, spin_json)[m
 [m
[31m-    print("g-matrix: ")[m
[31m-    print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\[m
[31m-                for row in np.round(np.sqrt(gmatrix[:,:]), 6)]))[m
[31m-    print()[m
[32m+[m[32m    # print("g-matrix: ")[m
[32m+[m[32m    # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\[m
[32m+[m[32m    #             for row in np.round(np.sqrt(gmatrix[:,:]), 6)]))[m
[32m+[m[32m    # print()[m
 [m
 if excitanalysis_gvalue_cut != 0:[m
     get_gtensor_analysis(energies_json, excitenergies_json, spin_json, soc_json, orbitmoment_json, [m
