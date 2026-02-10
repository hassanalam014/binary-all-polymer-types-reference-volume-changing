
import os,sys,math,csv,matplotlib.pyplot as plt,numpy as npy
# from p_params import *
# from s_params import *
from loadExperimentalData import *
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
from loadPhysicalConstants import *
from calculatePureVariables import calculateNewMolecularParameters
# from wrapperFunctions import calculateBinarySolubilitySwelling
# from calculateBinaryResidual import calculateBinarySSQ
from checkResults import *
from sympy import *
import warnings
import cmath
from scipy.optimize import bisect,fsolve
from scipy.interpolate import interp1d
from Parameters_of_Different_Polymers import *
from All_Functions import *
from Parameters_for_Mixtures_and_Tg import *
from scipy.interpolate import interp1d
from scipy.interpolate import make_interp_spline, BSpline
from collections import OrderedDict			#For Exotic Line Styles

###########################################################################
###########################################################################

def discard_zeros(x,y):
	
	for i in range(len(x)):
		if x[i]==0:
			y[i]=0

	for i in range(len(y)):
		if y[i]==0:
			x[i]=0

	x = npy.delete(x, npy.argwhere( (x >= 0) & (x <= 0) ))
	y = npy.delete(y, npy.argwhere( (y >= 0) & (y <= 0) ))
				
	return x,y

Polymer_Type='PMMA'
Solvent='CO2'

Parameters_Paper ='Self_Grassia'			# P*T*R* and g,epsilon_2,x (PVT-Tg Data Paper or Direct P*T*R* Values Reference)
Cp_Polymer_Weight = '02kilo_POST_THESIS'	# g,epsilon_2,x (Cp Paper Reference)
Paper_Number = 'Paper15'						# Solubility or Swelling Data Reference
#for PS: 'Paper4_11_12'
#for PMMA: 'Paper15'
kwargs = {'Polymer_Type':Polymer_Type,'Solvent':Solvent,'Parameters_Paper':Parameters_Paper,'Paper_Number':Paper_Number,'Cp_Polymer_Weight':Cp_Polymer_Weight}

# Ppstar,Tpstar,Rpstar,Mp,Psstar,Tsstar,Rsstar,Ms,P_exp,Tg_exp=Parameters_of_Different_Polymers(**kwargs)

forward=False
backward=True

Condo_Original=False 
CondoFind_Tg_at_P=False
CondoFind_Pg_at_T=False

Kier=False
Hassan=True  
Condo=False 

Find_Tg_at_P=False
Find_Pg_at_T=True




if False: #To hide irrelevant lines:

	############################################## PRE-THESIS Data #############################################

	# For zeta=1.124 and delta=0.97
	# P=[0.101325,0.75673333,1.41214167,2.06755,2.72295833]
	# Tg_bisect_fwd=[372.0703125,362.6953125,252.1484375,270.5078125,292.7734375]
	# Pg_bisect_fwd=[1.33789062,1.92382812,2.48242188,2.84179688,2.84960938,2.42382812,1.58398438]
	# T=[250.0,266.66666667,283.33333333,300.0,316.66666667,333.33333333,350.0]

	# For zeta=1.1050 and delta=0.97
	# P = [0.101325 ,1.75673333 ,3.41214167 ,5.06755 ,6.72295833 ,8.37836667 ,10.033775 ,11.68918333 ,13.34459167 ,15.0]
	# Tg_bisect_fwd = [373.2421875 ,365.8203125 ,359.5703125 ,353.3203125 ,347.8515625 ,343.1640625 ,339.2578125 ,336.1328125 ,334.9609375 ,333.7890625]
	# Pg_bisect_fwd = [1.54882812]
	# T = [366.66666667]

	# For zeta=1.100 and delta=0.97
	# P1 = [2.84570312,3.54301184,4.11662632,4.69024079,5.26385526,5.83746974,6.41108421,6.98469868,7.55831316,8.13192763,8.70554211,9.27915658,9.85277105,10.42638553,11.0000000]								
	# Tg1 = [265.78947368,265.4296875,266.6015625,265.8203125,265.0390625,264.2578125,263.4765625,262.6953125,261.9140625,261.1328125,260.7421875,259.9609375,259.1796875,258.3984375,257.6171875]								
	# P = [0.101325,0.67493947,1.24855395,1.82216842,2.39578289,2.96939737]							
	# Tg_bisect_fwd = [372.8515625,366.6015625,361.1328125,355.6640625,349.8046875,343.5546875]					
	# Pg_bisect_fwd = [1.87304688,2.84570312,3.32421875,3.77734375,4.16796875,4.44140625,4.58203125,4.55859375,4.36328125,4.00390625,3.49609375,2.86523438,2.12695312,1.32226562]								
	# T = [250.000000,265.78947368,273.68421053,281.57894737,289.47368421,297.36842105,305.26315789,313.15789474,321.05263158,328.94736842,336.84210526,344.73684211,352.63157895,360.52631579]								

	# # For zeta=1.088 and delta=0.97 #Full Data
	# P1 = [0.101325,3.42621458,6.75110417,10.07599375,13.40088333,16.72577292,20.0506625,23.37555208,26.70044167,30.02533125,33.35022083,36.67511042,40.00000]
	# P2 = [0.101325,1.09211364,2.08290227,3.07369091,4.06447955,5.05526818,6.04605682,7.03684545,8.02763409,9.01842273,10.00921136,11.00000]	
	# P=P1+P2

	# Tg_bisect_fwd1 = [372.8515625,345.8984375,294.7265625,291.6015625,288.8671875,286.1328125,283.7890625,281.4453125,279.4921875,277.5390625,275.1953125,273.2421875,271.2890625]
	# Tg_bisect_fwd2 = [372.8515625,364.6484375,356.8359375,348.6328125,270.5078125,329.1015625,290.4296875,294.3359375,293.5546875,292.3828125,291.6015625,290.8203125]	
	# Tg_bisect_fwd=Tg_bisect_fwd1+Tg_bisect_fwd2

	# Pg_bisect_fwd1 = [1.87304688,3.95703125,5.91015625,5.34765625,2.93359375]
	# Pg_bisect_fwd2 = [2.79882812,5.53515625,5.72265625,3.46484375]		
	# Pg_bisect_fwd=Pg_bisect_fwd1+Pg_bisect_fwd2

	# T1 = [250.0,275.0,300.0,325.0,350.0]					
	# T2 = [263.63636364,290.90909091,318.18181818,345.45454545]									
	# T=T1+T2

	# For zeta=1.088 and delta=0.97 #Data Upto 11 MPa without retrograde
	# P = [0.101325,1.09211364,2.08290227,3.07369091,3.42621458,5.05526818,5.34765625,5.72265625,5.91015625,6.04605682,6.75110417,7.03684545,8.02763409,9.01842273,10.00921136,11.00000]	
	# Tg_bisect_fwd = [372.8515625,364.6484375,356.8359375,348.6328125,345.8984375,329.1015625,325.0,318.18181818,300.0,290.4296875,294.7265625,294.3359375,293.5546875,292.3828125,291.6015625,290.8203125]	

	# For zeta=1.088 and delta=0.97 #Data Upto 11 MPa with retrograde
	# P = [0.101325,1.09211364,2.08290227,3.07369091,3.42621458,5.05526818,5.34765625,5.72265625,5.91015625]
	# Tg_bisect_fwd = [372.8515625,364.6484375,356.8359375,348.6328125,345.8984375,329.1015625,325.0,318.18181818,300.0]

	# Pg_bisect_fwd = [1.87304688,2.79882812,3.95703125,5.53515625]
	# T = [250.0,263.63636364,275.0,290.90909091]

	# For PC CO2 Kier Method Fitted zeta=1.06672 and delta=1.0
	# P = [0.101325,2.25110417,6.5506625,8.70044167,10.85022083,13.00]
	# Tg_bisect_fwd = [423.046875,413.671875,398.046875,391.015625,384.9609375,379.1015625]
	# Pg_bisect_fwd = [5.96484375]
	# T = [400.0]

	# P = [0.101325,2.25110417,6.5506625,8.70044167,10.85022083,13.0]
	# Tg_bisect_fwd = [377.5390625,368.1640625,338.0859375,302.9296875,298.2421875,293.9453125]
	# Pg_bisect_fwd = [5.23828125]
	# T = [350.0]

	############################################## POST-THESIS Data ############################################

	'''
	#1 PMMA-CO2 for Data of Solubility Chul Park Paper 11, PVT Self_Grassia, Cp 01kilo_POST_THESIS, zeta=1.015, delta=1.0 
	Pg_bisect_bwd = [4.41796875,4.03515625,3.26953125,2.18945312]
	T_bwd = [296.15384615,307.69230769,319.23076923 ,330.76923077]
	Pg_bisect_fwd = [2.53320312,3.36328125,3.98828125,4.38671875,4.41796875,4.03515625,3.26953125,2.18945312]
	T_fwd = [250.0,261.53846154,273.07692308,284.61538462,296.15384615,307.69230769,319.23076923,330.76923077]
	Pg = Pg_bisect_fwd + Pg_bisect_bwd
	T = T_fwd + T_bwd

	#2 PMMA-CO2 for Data of Solubility Chul Park Paper 11, PVT Self_Grassia, Cp 02kilo_POST_THESIS, zeta=1.015, delta=1.0 
	Pg_bisect_fwd = [1.62695312,4.06640625,3.69921875,2.98828125,2.00585938]
	T_fwd = [284.61538462,296.15384615,307.69230769,319.23076923,330.76923077]
	Pg_bisect_bwd = [3.98828125,2.52148438]
	T_bwd = [300.0,325.0]
	Pg = Pg_bisect_fwd + Pg_bisect_bwd
	T = T_fwd + T_bwd

	#3 PMMA-CO2 for Data of Solubility Chul Park Paper 11, PVT Self_Wen, Cp 02kilo_POST_THESIS, zeta=1.04, delta=1.0 
	Pg_bisect_fwd = [5.34765625,6.07421875,5.18359375,3.62109375,1.68554688]
	T_fwd = [300.0,316.66666667,333.33333333,350.0,366.66666667]
	Pg_bisect_bwd = [5.34765625,6.07421875,5.18359375,3.62109375,1.68554688]
	T_bwd = [300.0,316.66666667,333.33333333,350.0,366.66666667]
	Pg = Pg_bisect_fwd + Pg_bisect_bwd
	T = T_fwd + T_bwd

	#4 PS-CO2 for Data of Solubility Paper 6, PVT Self_Grassia, Cp 02kilo_POST_THESIS, zeta=0.975, delta=1.0 
	P = [0.101325,1.41777083,2.73421667,4.0506625,5.36710833,6.68355417,8.0]
	Tg = [353.7109375,343.5546875,334.1796875,324.8046875,299.8046875,299.4140625,299.0234375]
	'''
	'''
	indep_list, dep_list = zip(*sorted(zip(indep_list, dep_list)))
	indep_list=list(indep_list)
	dep_list=list(dep_list)
	indep_list=indep_list[::-1] # end to beginning, counting down by 1
	dep_list=dep_list[::-1]     # end to beginning, counting down by 1
	'''
	# P1 = [1.87304688,2.79882812,3.95703125,5.53515625,5.91015625,5.72265625,6.75110417,7.03684545,8.02763409,9.01842273,10.00921136,11.00000]	
	# Tg1 = [250.0,263.63636364,275.0,290.90909091,300.0,318.18181818,294.7265625,294.3359375,293.5546875,292.3828125,291.6015625,290.8203125]	

	# indep_list = indep_list+Tg1
	# dep_list = dep_list+P1

	# plt.plot(P_exp,Tg_exp,color='b',marker='o',ls='',label='Tg_exp_condo',ms=markersize)

	pass

if Polymer_Type == 'PS':
	#WRONG: PS/CO2 Self_Grassia_02kilo_POST_THESIS_Paper4_11_12:
	# P = [0.101325,		0.36612127,		0.87980179,		0.91640181,		1.51449041,		1.65827857,		2.14672763,		2.43675536,		2.80238104,		3.21523214,		3.47169389,		3.99370893,		4.15856889,		4.77218571,		4.86886009,		5.5506625,		5.59476207,		6.32913929,		6.36749645,		7.10761607,		7.24583339,		7.69166667,		7.88609286,		8.66456964,		9.44304643,		10.22152321,	11]
	# Tg = [353.7109375,	351.42857143,	347.4609375,	347.14285714,	342.85714286,	341.9921875,	338.57142857,	336.5234375,	334.28571429,	331.4453125,	330,			326.7578125,	325.71428571,	322.0703125,	321.42857143,	317.3828125,	317.14285714,	313.0859375,	312.85714286,	309.1796875,	308.57142857,	306.42857143,	306.4453125,	306.4453125,	306.0546875,	306.0546875,	306.0546875]
	#WRONG: Removing close points:
	P = [0.101325,		0.87980179,		1.51449041,		2.14672763,			2.80238104,		3.47169389,		4.15856889,		4.86886009,		5.5506625,		6.32913929,		7.10761607,		7.88609286,		8.66456964,		9.44304643,		10.22152321,	11]
	Tg = [353.7109375,	347.4609375,	342.85714286,	338.57142857,		334.28571429,	330,			325.71428571,	321.42857143,	317.3828125,	313.0859375,	309.1796875,	306.4453125,	306.4453125,	306.0546875,	306.0546875,	306.0546875]
	
	#CORRECT:
	Pg = [11.000,10.009,9.018,9.000,8.727,8.455,8.275,8.182,8.028,7.909,7.636,7.500,7.429,7.400,7.364,7.357,7.353,7.281,7.106,6.906,6.932,6.944,6.944,6.947,6.956,6.956,6.968,6.974,6.981,6.980,6.991,7.003,7.006,7.021,7.031,7.044,7.065,7.068,7.081,7.106,7.119,7.116,7.144,7.163,7.181,7.191,7.211,7.219,7.244,7.256,7.258,7.281,7.319,7.306,7.305,7.294,7.286,7.256,7.244,7.228,7.214,7.181,7.143,7.119,7.091,7.071,7.056,7.037,7.022,7.000,6.994,6.931,6.834,6.818,6.666,6.545,6.478,6.309,6.273,6.141,6.055,6.046,6.000,5.551,5.107,5.055,4.803,4.064,3.623,3.074,2.826,2.479,2.083,1.932,1.409,1.092,0.424,0.101]
	T = [306.021,306.099,306.211,306.224,306.259,306.293,306.323,306.345,306.367,306.380,306.450,306.465,306.493,306.495,306.502,306.507,306.512,306.526,306.571,306.658,306.755,306.789,306.792,306.818,306.838,306.857,306.885,306.911,306.921,306.932,306.974,307.000,307.053,307.085,307.143,307.184,307.273,307.276,307.316,307.429,307.447,307.467,307.579,307.641,307.712,307.727,307.832,307.842,307.974,308.000,308.023,308.105,308.237,308.286,308.300,308.368,308.396,308.500,308.571,308.636,308.701,308.857,309.021,309.143,309.262,309.354,309.429,309.515,309.545,309.688,309.714,310.000,310.455,310.564,311.364,311.936,312.273,313.182,313.377,314.091,314.545,314.615,314.870,317.405,320.000,320.296,321.818,326.311,329.091,332.539,334.124,336.364,339.001,340.000,343.636,345.854,350.909,353.811]

if Polymer_Type == 'PMMA':
	# WRONG: #PMMA/CO2 Self_Grassia_02kilo_POST_THESIS_Paper15:
	# Pg = [2.05908203,	1.57226562,	0.51000977,	2.32226562,	2.32275391,	2.75585938,	4.13681641,	4.65234375,	4.91464844,	4.83984375,	4.61933594,	4.44140625,	4.08408203,	3.79296875,	3.35546875,	3.35107422,	2.94140625,	2.46513672,	1.95507812,	1.48032227,	0.5574707,	0.49418945,	0.37817383,	0.20744629]	
	# T = [272.5,	273.57142857,	280,	287.5,	287.50001,	289.28571429,	295,	297.14285714,	302.5,	305,	310,	312.85714286,	317.5,	320.71428571,	325,	325.0001,	328.57142857,	332.5,	336.42857143,	340,	347,	347.5,	348.44444444,	349.88888889]
	# WRONG: #Removing third line:
	# Pg = [0.51000977,	2.32226562,	2.32275391,	2.75585938,	4.13681641,	4.65234375,	4.91464844,	4.83984375,	4.61933594,	4.44140625,	4.08408203,	3.79296875,	3.35546875,	3.35107422,	2.94140625,	2.46513672,	1.95507812,	1.48032227,	0.20744629]	
	# T = [280,	287.5,	287.50001,	289.28571429,	295,	297.14285714,	302.5,	305,	310,	312.85714286,	317.5,	320.71428571,	325,	325.0001,	328.57142857,	332.5,	336.42857143,	340,	349.88888889]
	# WRONG: #Removing close points:
	Pg = [2.05908203,	1.57226562,	0.51000977,	2.32226562,	2.32275391,	2.75585938,	4.13681641,	4.65234375,	4.91464844,	4.83984375,	4.61933594,	4.44140625,	4.08408203,	3.79296875,	3.35546875,	3.35107422,	2.94140625,	2.46513672,	1.95507812,	1.48032227,	0.20744629]	
	T = [272.5,	273.57142857,	280,	287.5,	287.50001,	289.28571429,	295,	297.14285714,	302.5,	305,	310,	312.85714286,	317.5,	320.71428571,	325,	325.0001,	328.57142857,	332.5,	336.42857143,	340,	349.88888889]

	#CORRECT DOME
	Pg = [0.101,0.387,0.412,0.722,1.033,1.124,1.343,1.654,1.899,1.964,2.275,2.585,2.636,2.895,3.206,3.318,3.516,3.827,3.908,4.137,4.387,4.448,4.622,4.678,4.716,4.719,4.758,4.772,4.810,4.847,4.885,4.903,4.922,4.922,4.941,4.960,4.758,4.660,4.448,4.360,4.041,3.963,3.722,3.404,3.066,2.747,2.410,2.120,2.054,2.000,1.716,1.500,1.360,1.094,1.004,1.001,0.648,0.501,0.292,0.200,0.150,0.129,0.101,0.051,0.001]
	T = [350.869,348.421,348.154,345.713,343.350,342.632,341.025,338.701,336.842,336.357,333.975,331.533,331.053,328.975,326.318,325.263,323.467,320.361,319.474,316.885,313.684,312.764,310.000,308.966,307.931,307.895,307.178,306.897,305.862,304.828,303.793,302.759,302.105,301.724,300.690,299.655,298.936,298.621,297.881,297.586,296.552,296.316,295.517,294.483,293.448,292.414,291.379,290.526,290.345,290.166,289.310,288.682,288.276,287.500,287.241,287.217,286.207,285.791,285.172,284.937,284.796,284.737,284.657,284.517,284.380]
	#Correct Dome Lower Part:
	Pg_lower = [0.001,0.051,0.101,0.101,0.150,0.200,0.412,0.501,0.722,0.894,1.001,1.033,1.343,1.500,1.654,1.770,1.964,2.000,2.275,2.585,2.895,3.206,3.516,3.827,4.137]
	T_lower = [276.842,276.739,276.636,276.631,276.533,276.430,275.986,275.811,275.361,275.000,274.775,274.717,274.053,273.721,273.389,273.158,272.705,272.627,272.002,271.260,270.459,269.600,268.643,267.568,266.201]

	#CORRECT PMMA Different zeta (0.88,0.90,0.92,0.93) but fixed delta = 0.88621038 (as iterated in Paper15)
	P_88 = [0.101,0.885,1.670,1.757,2.237,2.454,3.238,3.412,4.022,4.806,5.068,5.310,5.590,6.374,6.723,7.159,7.943,8.378,8.727,9.511,9.745,10.034,10.295,11.079,11.689,11.863,12.648,13.345,13.432,14.216,15.000]
	T_88 = [351.444,348.490,345.976,345.718,344.286,343.688,341.540,341.078,339.509,337.597,336.985,336.429,335.803,334.117,333.418,332.570,331.152,330.443,329.906,328.842,328.571,328.262,328.015,327.478,327.274,327.231,327.070,326.973,326.952,326.866,326.780]
	P_90 = [0.101,0.622,1.143,1.664,2.185,2.706,3.227,3.748,4.269,4.790,5.311,5.832,6.353,6.874,7.395,7.916,8.000,8.437,8.778,8.958,9.185,9.479,9.556,10.000,10.333,11.111,11.889,12.667,13.444,14.222,15.000]
	T_90 = [351.299,348.701,346.436,344.326,342.295,340.322,338.408,336.533,334.697,332.900,331.123,329.404,327.705,326.084,324.521,323.037,322.803,321.670,320.850,320.459,320.000,319.482,319.365,319.014,318.896,318.643,318.447,318.291,318.135,317.998,317.881]
	P_92 = [0.101,0.522,0.808,1.515,1.571,2.222,2.688,2.930,3.637,3.805,4.344,4.887,5.051,5.758,5.901,6.465,6.602,6.680,6.742,6.805,6.812,6.867,6.914,6.977,7.000,7.039,7.102,7.164,7.172,7.227,7.289,7.336,7.398,7.461,7.445,7.336,7.227,7.225,7.133,7.070,7.742,7.879,8.000,8.586,8.602,9.000,9.293,9.602,10.000,10.758,11.000,12.000,12.055]
	T_92 = [351.123,348.421,346.768,342.939,342.632,339.248,336.842,335.615,331.963,331.053,328.213,325.263,324.346,320.283,319.474,315.928,315.000,314.583,314.167,313.750,313.684,313.333,312.917,312.500,312.383,312.083,311.667,311.250,311.201,310.833,310.417,310.000,309.583,309.167,308.750,308.333,307.917,307.895,307.500,307.083,306.667,306.592,306.539,306.260,306.250,306.070,305.967,305.833,305.687,305.417,305.336,305.023,305.000]
	P_93 = [0.101,0.372,1.093,2.084,2.123,3.076,3.874,4.068,5.059,5.330,5.928,5.984,6.000,6.041,6.051,6.097,6.153,6.191,6.247,6.284,6.291,6.284,6.097,5.909,5.703,5.667,5.497,5.333,5.309,5.103,5.069,5.028,5.028,4.983,4.953,4.908,4.863,4.818,4.773,4.743,4.698,4.743,4.878,5.043,5.084,5.234,5.403,5.459,5.572,5.741,5.928,6.116,6.303,6.333,6.491,6.667,7.000,7.042,7.333,7.667,7.709,8.000,8.034,9.025,10.017,11.008,12.000]
	T_93 = [351.025,349.000,344.322,338.274,338.000,332.205,327.000,325.738,318.369,316.000,310.000,309.292,309.170,308.583,308.529,307.875,307.167,306.458,305.750,305.042,305.000,304.333,303.625,302.917,302.208,302.078,301.500,300.909,300.792,300.083,300.000,299.875,299.857,299.714,299.571,299.429,299.286,299.143,299.000,298.857,298.714,298.571,298.429,298.286,298.250,298.125,298.000,297.958,297.875,297.750,297.625,297.500,297.375,297.350,297.250,297.150,296.951,296.927,296.752,296.566,296.542,296.393,296.369,295.864,295.380,294.929,294.500]


if Find_Tg_at_P:
	#For no retrograde indep_list is P and dep_list is Tg
	dep_list = P#+Tg_bisect_fwd
	indep_list = Tg#_bisect_fwd#+P

if Find_Pg_at_T:
	#For retrograde indep_list is T and dep_list is Pg
	indep_list = T#+Tg_bisect_fwd
	dep_list = Pg#_bisect_fwd#+P

print indep_list
print dep_list

#Setting font size
axis_size = 24  #Size  of x and y axis wordings (names)
title_size = 20 #We have no title
size = 18		#Size of legendre
label_size = 22	#Size of values of ticks on x and y axis

plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

#Markers
mark1 = 'o'
mark2 = 's'
mark3 = '>'
mark4 = '<'
mark5 = 'D'
mark6 = 'H'
mark7 = 'P'
mark8 = 'X'

mark9 = '^'
mark10 = '1'
mark11 = '2'
mark12 = '3'
mark13 = '*'
mark14 = 'h'
mark15 = '+'

linestyles = OrderedDict(
    [('solid',               (0, ())),
     ('loosely dotted',      (0, (1, 10))),
     ('dotted',              (0, (1, 5))),
     ('densely dotted',      (0, (1, 1))),

     ('loosely dashed',      (0, (5, 10))),
     ('dashed',              (0, (5, 5))),
     ('densely dashed',      (0, (5, 1))),

     ('loosely dashdotted',  (0, (3, 10, 1, 10))),
     ('dashdotted',          (0, (3, 5, 1, 5))),
     ('densely dashdotted',  (0, (3, 1, 1, 1))),

     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])

#Linestyles
ls1 = '-'
ls2 = '--'
ls3 = ':'
ls4 = '-.'
ls5 = linestyles['densely dashdotdotted']

#Setting saved image properties
img_extension = '.png'
img_dpi = None
output_folder = 'Post_Thesis'

#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)

#General line properties.
linewidth = 3
markersize = 10

#==================================================================================
#Plots.
figPUREPS=plt.figure(num=None, figsize=(10,6), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()

if Polymer_Type == 'PMMA':
	#PMMA-CO2 Data Below:
	#Experiment Data PMMA from Condo Paper;
	Tg_exp1=	[378.0,277.8,297.7,317.7,337.8,347.9]
	P_exp1=	[0.101325,3.75,5.14,5.70,5.09,3.73]
	#Experiment Data PMMA from Ruosong et al.
	Tg_exp2 = [368.3,352.3,346.9,326.5,323.2,311.4,317.9,320.1,320.9,322.2,323.1,296.2,315.6]
	P_exp2 = [2,4,6,8,10,12,14,16,18,20,22,4,6]
	#Experiment Data PMMA from Maria Pantoula's Paper Ref. Kamiya [24]
	Tg_exp3 = [358.138,348.031,338.091,308.019]
	P_exp3 = [3.60067,3.90142,4.89557,5.00418]
	#Experiment Data PMMA from Maria Pantoula's Paper Ref. Wissinger[18]
	Tg_exp4 = [331.826,314.618,305.68]
	P_exp4 = [3.90142,3.96825,3.90142]
	#Experiment Data PMMA from Maria Pantoula's Paper Ref. Alessi[34]
	Tg_exp5 = [312.697,307.434]
	P_exp5 = [7.96992,10.0167]
	#Experiment Data PMMA from Maria Pantoula's Paper Ref. O'Neill[36]
	Tg_exp6 = [364.905,352.291,345.692,341.181,330.907]
	P_exp6 = [0.85213,1.77945,   2.23893,2.92398,3.67586]
	#Experiment Data PMMA from Chul Parks's Paper Ref. Handa
	Tg_exp7 = [367.2229,358.135,349.1906,346.4349,341.2781,333.2876,325.3884,313.574,306.2017,298.0507,290.2496,275.85918,271.3099]
	P_exp7 = [0,1.516703528,2.793439058,3.283467023,4.14107169,4.898557125,5.401372305,5.77479546,5.413126005,4.894818233,4.148741993,3.295585493,2.796499073]
	#Experiment Data PMMA from Chul Park's Paper Ref. Condo
	Tg_exp8 = [349.2809,338.7512,330.6587,319.2883,308.4899,298.0621,288.2614,278.34445]
	P_exp8 = [3.655390568,5.140430033,5.917521855,5.872158653,5.96386791,5.147036423,4.548094215,3.652857443]
	#Experiment Data PMMA from Dehua Liu's Paper Ref. His Work
	Tg_exp9 = [388.597,383.74,380.03,374.769,359.1871,344.4148]
	P_exp9 = [0.0988235,1.50706,2.24412,3.01,4.36882,5.99941]
	#Experiment Data PMMA from Dehua Liu's Paper Ref. Handa
	Tg_exp10 = [363.8414,354.1281,342.3911,336.9274,332.6778,324.1787]
	P_exp10 = [0.0947059,0.916177,1.88588,2.31824,3.04706,3.80059]
	# Handa et al. Original Data. Paper Title: A New Technique for Measuring Retrograde Vitrification
	P_exp11 = [0,	1.52279316,	2.8602527625,	3.4549899825,	4.1991005175,	4.9488041925,	5.4080496225,	5.9240978475,	5.4111299025,	4.9507091025,	4.1929804875,	3.4519299675,	2.861519325]
	Tg_exp11 = [368.1392,	355.8907,	348.0299,	344.0221,	339.2131,	331.0007,	324.4673,	312.1602,	305.2321,	299.1078,	291.1441,	277.71224,	270.592]
	# Condo data taken from Handa et al. paper. Title: A New Technique for Measuring Retrograde Vitrification
	P_exp12 = [0,	3.748113075,	5.1685781175,	5.8560682425,	5.857618515,	5.9338453125,	5.17061475,	4.5553592175,	3.75202422]	
	Tg_exp12 = [378.204,	348.0411,	338.1163,	328.144,	317.972,	308.1833,	298.0398,	288.2806,	278.13667]		

	# PMMA-CO2 Experimental Data Plots Below:
	# plt.plot(P_exp1,Tg_exp1,color='b',marker=mark1,ls='',label='Condo',ms=markersize)
	# plt.plot(P_exp2,Tg_exp2,color='g',marker=mark2,ls='',label='Ruosong',ms=markersize)
	# plt.plot(P_exp3,Tg_exp3,color='r',marker=mark3,ls='',label='Kamiya [24] from Pantoula',ms=markersize)
	# plt.plot(P_exp4,Tg_exp4,color='c',marker=mark4,ls='',label='Wissinger[18] from Pantoula',ms=markersize)
	# plt.plot(P_exp5,Tg_exp5,color='m',marker=mark5,ls='',label='Alessi[34] from Pantoula',ms=markersize)
	# plt.plot(P_exp6,Tg_exp6,color='y',marker=mark6,ls='',label='ONeill[36] from Pantoula',ms=markersize)
	# plt.plot(P_exp7,Tg_exp7,color='k',marker=mark1,ls='',label='Handa form C.Park',ms=markersize)	#Taken form: C.Park's Paper
	# plt.plot(P_exp8,Tg_exp8,color='k',marker=mark8,ls='',label='Condo from C.Park',ms=markersize)
	# plt.plot(P_exp9,Tg_exp9,color='b',marker=mark9,ls='',label='Dehua Liu',ms=markersize)
	# I do NOT use this handa ==> plt.plot(P_exp10,Tg_exp10,color='k',marker=mark10,ls='',label='Handa from Dehua Liu',ms=markersize)
	plt.plot(P_exp11,Tg_exp11,color='k',marker=mark1,ls='',label='Handa',ms=markersize)	#Taken form: C.Park's Paper
	# plt.plot(P_exp12,Tg_exp12,color='k',marker=mark8,ls='',label='Condo',ms=markersize)

if Polymer_Type == 'PS':
	#PS-CO2 Data Below:
	#Experiment Data PS;
	Tg_exp=[373.0,328.3,308.2,305.0,304.5]	
	P_exp=[0.101325,5.90,7.04,8.36,10.38]
	# Condo et al.			
	Tg_exp_condo = [373.18,328.40,308.21,305.84,304.49]
	P_exp_condo = [0.100,6.072,7.139,8.419,10.440]
	# Pham el al. 90 nm				
	Tg_exp_pham90nm = [383.11,348.03,322.87,308.10,298.06]
	P_exp_pham90nm = [0.089,3.990,4.792,5.209,4.884]
	# Pham et al. 17 nm				
	Tg_exp_pham17nm = [362.47,347.58,322.65,307.98,297.61]
	P_exp_pham17nm = [0.089,1.278,3.218,4.071,3.675]
	# Wissinger Bulk				
	Tg_exp_wissinger = [373.86,338.67,324.00,308.89]	
	P_exp_wissinger = [0.100,3.645,4.864,6.072]	
	#Condo Data Unknown Source:
	Tg_exp6=[374.1,339.2,323.9,308.8]
	P_exp6=[0.008250,3.64,4.88,6.06]
	#Experiment Data PS from Maria Pantoula's Paper Ref. Zhang[32]
	Tg_exp7=	[375.498,372.188,367.621,360.073,354.183,346.24,337.037,333.065]
	P_exp7 = [0.362753,0.730141,1.4465,2.0206,3.11808,4.20183,4.8678,5.38213]
	#Experiment Data PS from Maria Pantoula's Paper Ref. Tsivintzelis[33]
	Tg_exp8=	[348.288,339.218,318.033]
	P_exp8 = [2.99885,4.00457,6.00232]
	#Experiment Data PS from Maria Pantoula's Paper Ref. O'Neill[36]
	Tg_exp9 =	[372.121,367.42,362.852,354.444,344.78,330.614,321.212]
	P_exp9 = [0.541896,0.909318,1.36397,2.27326,3.59581,5.26281,6.07571]
	#Experiment Data PS from Maria Pantoula's Paper Ref. Wissinger[18]
	Tg_exp10=	[377.218,338.289,323.326,308.297]
	P_exp10 = [0.0137675,3.65106,4.85895,6.0852]

	# PS-CO2 Experimental Data Plots Below:
	plt.plot(P_exp_condo,Tg_exp_condo,color='r',marker=mark1,ls='',label='Condo',ms=markersize)
	# plt.plot(P_exp_pham90nm,Tg_exp_pham90nm,color='r',marker='P',ls='',label='Pham',ms=markersize)
	# plt.plot(P_exp_pham17nm,Tg_exp_pham17nm,color='b',marker='^',ls='',label='Pham et al. - h~17nm',ms=markersize)
	plt.plot(P_exp_wissinger,Tg_exp_wissinger,color='g',marker=mark2,ls='',label='Wissinger',ms=markersize)
	# plt.plot(P_exp6,Tg_exp6,color='y',marker='2',ls='',label='Condo-2',ms=markersize)
	# plt.plot(P_exp7,Tg_exp7,color='b',marker='3',ls='',label='Zhang[32] from Pantoula',ms=markersize)
	# plt.plot(P_exp8,Tg_exp8,color='k',marker='*',ls='',label='Tsivintzelis[33] from Pantoula',ms=markersize)
	# plt.plot(P_exp9,Tg_exp9,color='b',marker='h',ls='',label='ONeill[36] from Pantoula',ms=markersize)
	# plt.plot(P_exp10,Tg_exp10,color='k',marker='+',ls='',label='Wissinger[18] from Pantoula',ms=markersize)

##################################
#CO2 Saturation Line:
##################################
#OHIO UNIVERSITY:
P_CO2_vap_ohio = [2,	2.1,	2.2,	2.3,	2.4,	2.5,	2.6,	2.7,	2.8,	2.9,	3,	3.1,	3.2,	3.3,	3.4,	3.5,	3.6,	3.7,	3.8,	3.9,	4,	4.1,	4.2,	4.3,	4.4,	4.5,	4.6,	4.7,	4.8,	4.9,	5,	5.1,	5.2,	5.3,	5.4,	5.5,	5.6,	5.7,	5.8,	5.9,	6,	6.1,	6.2,	6.3,	6.4,	6.5,	6.6,	6.7,	6.8,	6.9,	7,	7.1,	7.2,	7.3,	7.377]
T_CO2_vap_ohio = [253.65,	255.25,	256.79,	258.29,	259.73,	261.14,	262.5,	263.83,	265.12,	266.37,	267.6,	268.79,	269.96,	271.1,	272.22,	273.31,	274.38,	275.43,	276.45,	277.46,	278.45,	279.42,	280.37,	281.31,	282.23,	283.13,	284.02,	284.89,	285.75,	286.6,	287.43,	288.26,	289.06,	289.86,	290.65,	291.42,	292.18,	292.93,	293.68,	294.41,	295.13,	295.84,	296.54,	297.23,	297.92,	298.59,	299.26,	299.92,	300.56,	301.2,	301.83,	302.45,	303.07,	303.67,	304.13]
plt.plot(P_CO2_vap_ohio,T_CO2_vap_ohio,color='k',marker='',ms=markersize,lw=1.0,ls='--')#,label='Saturation ')
#plt.plot(7.377,304.13,color='k',marker='o',ms=4,lw=2,ls='')	#Critical point

#CONDO 1994:
# P_CO2_vap_condo1994 =	[3.43058,		3.57993,	3.745,		3.96518,	4.15403,	4.31134,	4.47644,	4.71247,	4.89342,	5.07426,	5.23955,	5.44413,	5.62515,	5.83762,	6.09747,	6.31784,	6.51457,	6.72709,	6.96325,	7.21517,	7.42769,	7.65596,	7.89216,	8.11256,	8.3882,	8.59285,	8.76595,	8.88393,	8.98628]
# T_CO2_vap_condo1994 =	[273.0575061,	274.94582,	277.02298,	279.41536,	281.0529,	282.62704,	284.5785,	286.7197,	288.3571,	290.4345,	291.6945,	293.458,	294.844,	296.6076,	298.1206,	299.8216,	301.4592,	303.0343,	304.6726,	306.3741,	307.9491,	309.5873,	311.0999,	312.6751,	314.1256,	315.6376,	317.1492,	318.3455,	319.0387]
#Condo Smooth out:
# P_CO2_vap_condo1994 =	[3.43,		3.71,	3.99,	4.26,	4.54,	4.82,	5.1,	5.38,	5.65,	5.93,	6.21,	6.49,	6.76,	7.04,	7.32,	7.6,	7.88,	8.15,	8.43,	8.71,	8.99]
# T_CO2_vap_condo1994 =	[273.89,	276.72,	279.48,	282.18,	284.82,	287.39,	289.89,	292.33,	294.71,	297.02,	299.27,	301.45,	303.57,	305.62,	307.61,	309.53,	311.39,	313.18,	314.91,	316.58,	318.18]
#Condo Smooth out only above 300 K:
# P_CO2_vap_condo1994 =	[5.93,	6.21,	6.49,	6.76,	7.04,	7.32,	7.6,	7.88,	8.15,	8.43,	8.71,	8.99]
# T_CO2_vap_condo1994 =	[297.02,	299.27,	301.45,	303.57,	305.62,	307.61,	309.53,	311.39,	313.18,	314.91,	316.58,	318.18]
# plt.plot(P_CO2_vap_condo1994,T_CO2_vap_condo1994,color='k',marker='',ms=markersize,lw=1,ls='-',label='Saturation ')
# #plt.plot(7.32,307.61,color='k',marker='o',ms=4,lw=2,ls='')		#Critical point

# plt.text(7.6,306,'Critical Point',horizontalalignment='right')	#Write at Critical point
##################################



#For both: No retrograde indep_list is P and dep_list is Tg as well as for retrograde indep_list is T and dep_list is Pg
plt.plot(dep_list,indep_list,color='r',marker='',ms=markersize,lw=linewidth,ls='-',label=r'$\zeta_{sp}$ = 0.944')

if Polymer_Type == 'PMMA':
	plt.plot(Pg_lower,T_lower,color='r',marker='',ms=markersize,lw=linewidth,ls='-')

'''
plt.plot(P_93,T_93,color='b',marker='',ms=markersize,lw=linewidth,ls=ls2,label='0.93')
plt.plot(P_92,T_92,color='g',marker='',ms=markersize,lw=linewidth,ls=ls3,label='0.92')
plt.plot(P_90,T_90,color='m',marker='',ms=markersize,lw=linewidth,ls=ls4,label='0.90')
plt.plot(P_88,T_88,color='k',marker='',ms=markersize,lw=linewidth,ls=ls5,label='0.88')
plt.axis([-0.3,17.5,265,355])
'''

if False:		#Old sorting plots
		
	# if Find_Tg_at_P:
	# 	if forward:
	# 		P,Tg_bisect_fwd=discard_zeros(P,Tg_bisect_fwd)
	# 		plt.plot(P,Tg_bisect_fwd,color='k',marker='x',lw=linewidth,ls='-.',label='Tg_bisect_fwd')
	# 	if backward:
	# 		P,Tg_bisect_bwd=discard_zeros(P,Tg_bisect_bwd)
	# 		plt.plot(P,Tg_bisect_bwd,color='b',lw=linewidth,ls='-',label='Tg_bisect_bwd')

	# if Find_Pg_at_T:
	# 	if forward:
	# 		Pg_bisect_fwd,T=discard_zeros(Pg_bisect_fwd,T)
	# 		plt.plot(Pg_bisect_fwd,T,color='k',marker='x',lw=linewidth,ls='-.',label='Pg_bisect_fwd')
	# 	if backward:
	# 		Pg_bisect_bwd,T=discard_zeros(Pg_bisect_bwd,T)
	# 		plt.plot(Pg_bisect_bwd,T,color='b',lw=linewidth,ls='-',label='Pg_bisect_bwd')
	pass

plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel(r'Glass Temperature Tg (K)',fontsize=axis_size)
# plt.axis([-0.3,11.5,290,390])
plt.axis([-0.3,6.5,265,370])
plt.legend(loc=1,fontsize=size,numpoints=1,frameon=True)
# plt.title(kwargs, fontdict=None, loc='center', pad=None)
plt.subplots_adjust(left=0.15,right=0.95,top=0.97,bottom=0.12,wspace=0.30,hspace=0.25)
figPUREPS.savefig('./'+output_folder+r'\PMMA_CO2_Self_Grassia_02kilo_POST_THESIS_Paper15_Tg(P)_corr'+img_extension,dpi=240)

plt.show()
