
bd = '/mnt/c/scratch/sciteam/dcollins/Paper33_Galaxy'
weight_fields = {'scaled_div_b':'cell_volume','magnetic_field_strength':'density'}
methods = {'abs_divb':'mip'}
sim_base_dir = {}
centers = {}
sim_base_dir['e02']= 'e02_64_stock_nograckle'
sim_base_dir['e04']= 'e04_64_mhdct'  
sim_base_dir['h01']= 'h01_64_stock_newCrosby'
sim_base_dir['e06']='e06_crosby_mhdct'
sim_base_dir['e07']='e07_crosby_ppm'
sim_base_dir['e08']='e08_crosby_mhdct_bsmall'
sim_base_dir['e08b']='e08b_no_diskgrav_no_isolated'
sim_base_dir['e08c']='e08c_accelfix'
sim_base_dir['e08d']='e08d_enzo-dev'
sim_base_dir['e08e']='e08e_no_grackle_enzo-dev'
sim_base_dir['e08f']='e08f_unigrid'
sim_base_dir['e08g']='e08g_onelev'
sim_base_dir['e08h']='e08h_twolev'
sim_base_dir['e08i']='e08i_128root'
sim_base_dir['e08j']='e08j_wtf'
sim_base_dir['e08k']='e08k_e08i_long'
sim_base_dir['e09'] = 'e09_crosby_mhdct_bsmall_dev'
sim_base_dir['e10'] = 'e10_crosby_mhdct_bsmall_nograckle'
sim_base_dir['h02'] = 'h02_64_stock_crosby_grackle'
sim_base_dir['d02'] = 'd02_64_stock_crosby_dedner'
sim_base_dir['d03'] = 'd03_64_hydro3'
sim_base_dir['d04'] = 'd04_ppm_repeat'
sim_base_dir['i01'] = 'i01_64_crosby_mhd'
sim_base_dir['i02'] = 'i02_monkey_around'
sim_base_dir['i03'] = 'i03_monkey_around_ppm'
sim_base_dir['i04'] = 'i04_monkey_around_newPF_ppm'
sim_base_dir['i05'] = 'i05_i02_nofield'
sim_base_dir['a01'] = 'a01_agora1'
sim_base_dir['ai01'] = 'ai01_newDevin'
sim_base_dir['ai02'] = 'ai02_weak_field'
sim_base_dir['ai03'] = 'ai03_repeat_i02'
sim_base_dir['ai14'] = 'ai14_fixed_zerofield'
sim_base_dir['ai16'] = 'ai16_fixed_tinyfield'
sim_base_dir['ai17'] = 'ai17_repeat_i02'
sim_base_dir['ai17b'] = 'ai17b_rerun'

sim_base_dir['ai16b'] = 'ai16b'
sim_base_dir['ai21b'] = 'ai21b'
sim_base_dir['ai21b'] = 'ai21b'
sim_base_dir['ai22b'] = 'ai22b'
sim_base_dir['ai23b'] = 'ai23b'
sim_base_dir['ai24b'] = 'ai24b_longerdt'
centers['ai17'] = nar([ 0.49853515625, 0.49853515625, 0.49853515625])
centers['ai17b'] = nar([ 0.49853515625, 0.49853515625, 0.49853515625])
sim_base_dir['ai21'] = 'ai21_ai17_mass'
sim_base_dir['ai22'] = 'ai22_ai17_32cube'
sim_base_dir['ai23'] = 'ai23_length_units'
sim_base_dir['ai24'] = 'ai24_length_64'
sim_base_dir['ai25'] = 'ai25_smallroot_refine_on_start'
all_sim_list = ['ai16b','ai17b','ai21b','ai22b','ai23b','ai24b','i02','ai25']
all_fields_simple = [ "Cooling_Time", "Metal_Density", "GasEnergy", "Temperature", "TotalEnergy"]
all_fields_1 = ["Cooling_Time",          
"Electron_Density",      
"GasEnergy",             
"H2II_Density",          
"H2I_Density",           
"HII_Density",           
"HI_Density",            
"HM_Density",            
"HeIII_Density",         
"HeII_Density",          
"HeI_Density",           
"Metal_Density",         
"Temperature",           
"TotalEnergy"]
other_fields = ["Density",               
"x-velocity",            
"y-velocity",            
"z-velocity"]            
