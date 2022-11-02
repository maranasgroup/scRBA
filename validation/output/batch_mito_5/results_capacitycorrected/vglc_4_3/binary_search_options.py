path_gams = '../../../../../GAMS/'
path_pycore = '../../../../../pycore/'
path_enz_mw = '../../../../../input/enz_mw_g_per_mmol.txt'
path_pro_mw = '../../../../../input/pro_mw_g_per_mmol.txt'

report_file = './binary_search_report.txt' # Text file recording binary search process
mu_tol = 1e-5; # Tolerance of upper and lower bound gap to tolerance search
maxiter = 100; # Maximum number of iteration
mu_min0 = 0; mu_max0 = 0.5; # User-set initial upper and lower bounds of mu
biom_id = 'BIOSYN-BIODILAERO' # Set the ID of the biomass reaction
