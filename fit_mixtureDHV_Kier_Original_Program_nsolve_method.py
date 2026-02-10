# Date: May 2017
#
# Description	: The purpose of this file is to estimate the multicomponent fluid parameters
#				  for the PS/CO2 binary mixture.
#

import os,sys,math,csv,numpy as npy
# from p_params import *
# from s_params import *
from loadExperimentalData import *
from lmfit import minimize, Parameters, report_fit
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
# from calculateBinaryResidual import calculateBinaryResidualCHV
from Parameters_of_Different_Polymers import *
from SplittingExperimentalDataOf_X_Sw_aboveANDbelowTgANDIsotherms import *
from Parameters_of_Different_Polymers import *
from Split_Exp_Data_in_Isotherms import*
from loadPhysicalConstants import *
from checkResults import *
from sympy import *
import warnings

def binaryPhaseEquilibriumCHV(P,T,Mp,Ms,**kwargs):

	#Reference:
	# -p --> polymer
	# -s --> solvent

	for key,value in kwargs.iteritems():
		exec('{} = value'.format(key))
	
	if 'alpha_p' in kwargs and 'vhp' in kwargs and 'epsilon_p' in kwargs:
		Ppstar,Tpstar,Vpstar = calculateCharacteristicParameters(alpha_p,vhp,epsilon_p,Mp)
	elif 'Ppstar' in kwargs and 'Tpstar' in kwargs and 'Rpstar' in kwargs:
		pass
	else:
		raise ValueError('In binaryPhaseEquilibriumCHV, polymer parameters: Either molecular (alpha_p,vhp,epsilon_p) or characteristic (Ppstar,Tpstar,Rpstar) parameters must be passed into keyword arguments.')
	
	if 'alpha_s' in kwargs and 'vhs' in kwargs and 'epsilon_s' in kwargs:
		Psstar,Tsstar,Vsstar = calculateCharacteristicParameters(alpha_s,vhs,epsilon_s,Ms)
	elif 'Psstar' in kwargs and 'Tsstar' in kwargs and 'Rsstar' in kwargs:
		pass
	else:
		raise ValueError('In binaryPhaseEquilibriumCHV, solvent parameters: Either molecular (alpha_s,vhs,epsilon_s) or characteristic (Psstar,Tsstar,Rsstar) parameters must be passed into keyword arguments.')
	
	if 'k12' in kwargs and 'delta' in kwargs:
		pass
	elif 'zeta' in kwargs and 'delta' in kwargs:
		pass
	else:
		raise ValueError('In binaryPhaseEquilibriumCHV, mixture parameters: (k12,delta) or (zeta,delta) mixture parameters must be passed into keyword arguments.')
	
	#Allows for method argument in kwargs. Options are: 'disparate', 'single', 'mixed'.
	#Default option is 'disparate'.
	# -'disparate'	--> Mixture phase has constant hole volume of vhm. Pure phases have constant hole volumes vhp, vhs.
	# -'single'		--> All phases have hole volume vhm.
	method = kwargs.pop('method','disparate')
	
	#Boolean determining whether information is printed.
	#Default option is False.
	verbose = kwargs.pop('verbose',False)
	if verbose:
		print('FOR: P = {}MPa, T = {}K;'.format(P,T))
	
	#Initializing phi_p, phi_s and v_h as symbolic variables for the sympy package.
	#This is a step necessary for the numerical solver nsolve.
	phi_p = Symbol('phi_p',real=True)
	phi_s = Symbol('phi_s',real=True)
	v_h = Symbol('v_h',real=True)
	
	#PURE FLUID PARAMETERS.
	vpp = Mp/Rpstar				#Hassan: This is volume of whole polymer chain N_k*v_k. Not only segment of chain v_k.
	vss = Ms/Rsstar
	# vhp = kB*Tpstar/Ppstar
	# vhs = kB*Tsstar/Psstar
	vhm = delta*vhs											#Hassan: This is definition of delta
	# vhm=(phi_s*vhs+phi_p*vhp)/(phi_s+phi_p)

	# print 'vhp is:', vhp
	# print 'vhs is:', vhs

	chi_pp = -2*Tpstar/T			#Hassan: He has introduced a new variable not mentioned in paper.
	chi_ss = -2*Tsstar/T			#Hassan: He has introduced a new variable not mentioned in paper.
	alpha_p0 = (Ppstar*Mp)/(kB*Tpstar*Rpstar)
	alpha_s0 = (Psstar*Ms)/(kB*Tsstar*Rsstar)

	# vr = vhs
	alpha_p = alpha_p0*vhp/vr								#Hassan: v_h is volume of hole in mixture or pure solvent
	alpha_s = alpha_s0*vhs/vr
	alpha_0 = v_h/vr

	# print 'alpha_s0 is:', alpha_s0
	# print 'alpha_p0 is:', alpha_p0
	
	# vhs=vrs
	# Pfalse= 122.838243
	# vhp=kB*Tpstar/Pfalse

	#MIXTURE PARAMETERS.
	if 'k12' in kwargs:
		chi_ps = -(1.0-k12)*math.sqrt(chi_pp*chi_ss)
	elif 'zeta' in kwargs:
		chi_ps = -zeta*math.sqrt(chi_pp*chi_ss)				#Hassan: This is Kier paper Eq(34)
	
	# alpha_0 = v_h/vhs			#Hassan: Self Added Line
	# alpha_p0 = alpha_p0*vhp/vhs #Hassan: Self Added Line
	# alpha_p = alpha_p0			#Hassan: Self Added Line
	# alpha_s = alpha_s0			#Hassan: Self Added Line
	# vrp = kB*Tpstar/Ppstar
	# vrs = kB*Tsstar/Psstar
	# vhs=alpha_0*vrs
	# print 'vhs is:', vhs
	# print 'vrs is:', vrs
	# vhp=alpha_0*vrp

	# print 'am i running'
	#EQUATION OF STATE, CHEMICAL POTENTIAL.
	#Mixture equation of state.
	# EOS = v_h*P/(kB*T)-(chi_pp/2)*phi_p**2-chi_ps*phi_p*phi_s-(chi_ss/2)*phi_s**2+(1.0-1.0/alpha_p)*phi_p+(1.0-1.0/alpha_s)*phi_s+log(1.0-phi_p-phi_s)
	# EOS = v_h*P/(kB*T)-(chi_pp/2)*phi_p**2-chi_ps*phi_p*phi_s-(chi_ss/2)*phi_s**2+phi_p+(1.0-1.0/alpha_s)*phi_s+log(1.0-phi_p-phi_s)
	# EOS = vhs*P/(kB*T)-(chi_pp/2)*phi_p**2-chi_ps*phi_p*phi_s-(chi_ss/2)*phi_s**2+((1.0/alpha_0)-(1.0/alpha_p))*phi_p+((1.0/alpha_0)-(1.0/alpha_s))*phi_s+(1/alpha_0)*log(1.0-phi_p-phi_s)  #Hassan: My Modification of reference volume
	# EOS = vr*P/(kB*T)-(chi_pp/2)*phi_p**2-chi_ps*phi_p*phi_s-(chi_ss/2)*phi_s**2+((1.0/alpha_0)-(1.0/alpha_p))*phi_p+((1.0/alpha_0)-(1.0/alpha_s))*phi_s+(1/alpha_0)*log(1.0-phi_p-phi_s)  #Hassan: My Modification of reference volume
	EOS = vr*P/(kB*T)-(chi_pp/2)*phi_p**2-chi_ps*phi_p*phi_s-(chi_ss/2)*phi_s**2+((1.0/(vhp/vr))-(1.0/alpha_p))*phi_p+((1.0/(vhs/vr))-(1.0/alpha_s))*phi_s+(1/alpha_0)*log(1.0-phi_p-phi_s)  #Hassan: My Approximate Equation

	#Mixture equation of state in general.
	EOS_m = EOS.subs([(phi_p,phi_p),(phi_s,phi_s),(v_h,vhm)])
	
	#Mixture solvent chemical potential.
	# mu_s = alpha_s*(chi_ss*phi_s+chi_ps*phi_p+(1.0/alpha_s)*(1+log(phi_s))-log(1-phi_p-phi_s)-1)  #Hassan: This is my correction.
	# mu_s = (chi_ss*phi_s+chi_ps*phi_p+(1.0/alpha_s)*(1+log(phi_s))-(1/alpha_0)*log(1-phi_p-phi_s)-(1/alpha_0))  #Hassan: This is my correction when reference volume is not equal hole volume
	mu_s = alpha_s*(chi_ss*phi_s+chi_ps*phi_p+(1.0/alpha_s)*(1+log(phi_s))-(1/alpha_0)*log(1-phi_p-phi_s)-(1/alpha_0))  #Hassan: This is my correction when reference volume is not equal hole volume
	# mu_s = alpha_s*(chi_ss*phi_s+chi_ps*phi_p+(1.0/alpha_s)*(1+log(phi_s))-(1/alpha_0)*log(1-phi_p-phi_s)-(1/alpha_0)+(((vr*vhs)/(vhm**2))*((phi_p*(1-phi_p-phi_s))/((phi_p+phi_s)**2))))  #Hassan: vo derivative term added

	#Mixture solvent chemical potential in general.
	mu_s_m = mu_s.subs([(phi_p,phi_p),(phi_s,phi_s),(v_h,vhm)])

	#Mixture polymer chemical potential.
	# mu_p = alpha_p*(chi_pp*phi_p+chi_ps*phi_s+(1.0/alpha_p)*(1+log(phi_p))-log(1-phi_p-phi_s)-1) #Hassan: This is my correction
	# mu_p = (chi_pp*phi_p+chi_ps*phi_s+(1.0/alpha_p)*(1+log(phi_p))-(1/alpha_0)*log(1-phi_p-phi_s)-(1/alpha_0)) #Hassan: This is my correction when reference volume is not equal hole volume
	# mu_p = alpha_p*(chi_pp*phi_p+chi_ps*phi_s+(1.0/alpha_p)*(1+log(phi_p))-(1/alpha_0)*log(1-phi_p-phi_s)-(1/alpha_0)) #Hassan: This is my correction when reference volume is not equal hole volume

	#Mixture polymer chemical potential in general.
	# mu_p_m = mu_p.subs([(phi_p,phi_p),(phi_s,phi_s),(v_h,vhm)])
	
	#Selecting appropriate approach to use for limiting hole volumes.
	#Refers to 'method' in kwargs.
	if method == 'disparate':
		#Mixture equation of state in phi_s --> 0 limit. v_h takes on value vhp.
		EOS_p0 = EOS.subs([(phi_p,phi_p),(phi_s,0.0),(v_h,vhp)])
		#Mixture equation of state in phi_p --> 0 limit. v_h takes on value vhs.
		EOS_s0 = EOS.subs([(phi_p,0.0),(phi_s,phi_s),(v_h,vhs)])
		#Mixture solvent chemical potential in phi_p --> 0 limit. v_h takes on value vhs.
		mu_s0 = mu_s.subs([(phi_p,0.0),(phi_s,phi_s),(v_h,vhs)])
	if method == 'single':
		# print 'is this running'
		#Mixture equation of state in phi_s --> 0 limit. v_h takes on value vhm.
		EOS_p0 = EOS.subs([(phi_p,phi_p),(phi_s,0.0),(v_h,vhm)])
		#Mixture equation of state in phi_p --> 0 limit. v_h takes on value vhm.
		EOS_s0 = EOS.subs([(phi_p,0.0),(phi_s,phi_s),(v_h,vhm)])
		#Mixture solvent chemical potential in phi_p --> 0 limit. v_h takes on value vhm.
		mu_s0 = mu_s.subs([(phi_p,0.0),(phi_s,phi_s),(v_h,vhm)])
	if method == 'mixed':
		#Mixture equation of state in phi_s --> 0 limit. v_h takes on value vhp.
		EOS_p0 = EOS.subs([(phi_p,phi_p),(phi_s,0.0),(v_h,vhp)])
		#Mixture equation of state in phi_p --> 0 limit. v_h takes on value vhs.
		EOS_s0 = EOS.subs([(phi_p,0.0),(phi_s,phi_s),(v_h,vhs)])
		#Mixture solvent chemical potential in phi_p --> 0 limit. v_h takes on value vhm.
		mu_s0 = mu_s.subs([(phi_p,0.0),(phi_s,phi_s),(v_h,vhm)])

	print('P = {}, T = {};'.format(P,T))

	#CALCULATION OF PURE FLUID STATE AT P, T.
	#Solving for the volume fraction of the system in the pure fluid limiting cases.
	phip0 = nsolve(EOS_p0,phi_p,0.67,verify=True)
	phis0 = nsolve(EOS_s0,phi_s,0.1,verify=True)

	# print 'Is phip0 complex:',phip0
	# print 'Is phis0 complex:',phis0

	phip0=abs(phip0)
	phis0=abs(phis0)

	#CHECKING IF PURE VOLUME FRACTION RESULTS ARE VALID.
	checkVolumeFraction(phip0,'phi_p')
	checkVolumeFraction(phis0,'phi_s')
	
	#PRINTING OF RESULTS OF PURE FLUID CALCULATIONS.
	#FOR DIAGNOSTIC PURPOSES.
	if verbose:
		print('phip0 = {}, phis0 = {};'.format(phip0,phis0))
	
	# mu_gasPhase=mu_s0.subs(phi_s,phis0)
	# print mu_gasPhase


	#CALCULATION OF BINARY MIXTURE COMPOSITION AT P, T.
	#default [0.75,0.05]
	#Other good range [0.85,0.05]
	phip,phis = nsolve([EOS_m,(mu_s_m-mu_s0.subs(phi_s,phis0))],[phi_p,phi_s],[0.85,0.05],verify=True)#, prec=7)

	# if P<0.1:
	# 	phip,phis = nsolve([EOS_m,(mu_s_m-mu_s0.subs(phi_s,phis0))],[phi_p,phi_s],[0.90,0.000000002],verify=True)
	# elif P>0.7:
	# 	phip,phis = nsolve([EOS_m,(mu_s_m-mu_s0.subs(phi_s,phis0))],[phi_p,phi_s],[0.85,0.05],verify=True)
	
	# >>> nsolve(f, bounds(100), solver='bisect', verify=False)
	print 'Is phip complex:', phip
	print 'Is phis complex:', phis
	
	phip=abs(phip)
	phis=abs(phis)

	#CHECKING IF MIXTURE VOLUME FRACTION RESULTS ARE VALID.
	checkVolumeFraction([phip,phis],['phi_p','phi_s'])
	
	#PRINTING OF RESULTS OF MIXTURE COMPOSITION CALCULATIONS.
	#FOR DIAGNOSTIC PURPOSES.
	if verbose:
		print('phip = {}, phis = {};'.format(phip,phis))
	# phip0=0.0 #junk
	return [P,T,phip,phis,phip0]

def binarySolubilitySwellingCHV(P,T,Mp,Ms,**kwargs):
	# print 'this is also a great great great problem'
	for key,value in kwargs.iteritems():
		exec('{} = value'.format(key))
	
	if 'alpha_p' in kwargs and 'vhp' in kwargs and 'epsilon_p' in kwargs:
		Ppstar,Tpstar,Vpstar = calculateCharacteristicParameters(alpha_p,vhp,epsilon_p,Mp)
	elif 'Ppstar' in kwargs and 'Tpstar' in kwargs and 'Rpstar' in kwargs:
		pass
	else:
		raise ValueError('In binarySolubilitySwellingCHV, polymer parameters: Either molecular (alpha_p,vhp,epsilon_p) or characteristic (Ppstar,Tpstar,Rpstar) parameters must be passed into keyword arguments.')
	
	if 'alpha_s' in kwargs and 'vhs' in kwargs and 'epsilon_s' in kwargs:
		Psstar,Tsstar,Vsstar = calculateCharacteristicParameters(alpha_s,vhs,epsilon_s,Ms)
	elif 'Psstar' in kwargs and 'Tsstar' in kwargs and 'Rsstar' in kwargs:
		pass
	else:
		raise ValueError('In binarySolubilitySwellingCHV, solvent parameters: Either molecular (alpha_s,vhs,epsilon_s) or characteristic (Psstar,Tsstar,Rsstar) parameters must be passed into keyword arguments.')
	
	if 'k12' in kwargs and 'delta' in kwargs:
		pass
	elif 'zeta' in kwargs and 'delta' in kwargs:
		pass
	else:
		raise ValueError('In binarySolubilitySwellingCHV, mixture parameters: (k12,delta) or (zeta,delta) mixture parameters must be passed into keyword arguments.')
	
	#Boolean determining whether information is printed.
	#Default option is False.
	verbose = kwargs.get('verbose',False)
	
	# Boolean that determines method of calculation.
	#	 True: Uses simplified (original) swelling calculation assuming pure polymer.
	#	 False: Uses more sophisticated swelling calculation assuming air (N2) content.
	simplified = kwargs.pop('simplified',True)
	
	#PURE FLUID PARAMETERS.
	# vhp = kB*Tpstar/Ppstar
	# vhs = kB*Tsstar/Psstar
	alpha_p0 = (Ppstar*Mp)/(kB*Tpstar*Rpstar)
	alpha_s0 = (Psstar*Ms)/(kB*Tsstar*Rsstar)
	
	#MIXTURE PARAMETERS.
	# vhm = delta*vhs

	alpha_p = alpha_p0#*vhp/vhm
	alpha_s = alpha_s0#*vhs/vhm
	
	# CALCULATION OF VOLUME FRACTIONS AT P, T.
	if verbose:
		print('High-pressure solvent environment:')
	[Pd,Td,hsol_phip,hsol_phis,hphip0] = binaryPhaseEquilibriumCHV(P,T,Mp,Ms,**kwargs)
	
	#CALCULATION OF SOLVENT SOLUBILITY (MASS FRACTION) AT P, T.
	ms = (Ms*hsol_phis/alpha_s)/(Mp*hsol_phip/alpha_p+Ms*hsol_phis/alpha_s)   #Kier Original
	# ms = (Ms*hsol_phis/alpha_s0)/(Mp*hsol_phip/alpha_p0+Ms*hsol_phis/alpha_s0)   #Condo Solubility

	#CALCULATION OF VOLUME SWELLING (Vp/Vm) IN THE MIXTURE AT P, T.
	Sw = hphip0/hsol_phip

	#PRINTING OF RESULTS OF SOLUBILITY AND SWELLING.
	#FOR DIAGNOSTIC PURPOSES.
	if verbose:
		print('ms = {}, Sw = {};'.format(ms,Sw))
	
	return [P,T,ms,Sw]

def calculateBinarySolubilitySwelling(theory,P0,T0,Mp,Ms,**kwargs):

	if not isListOrNpyArray(P0) and not isListOrNpyArray(T0):
		exec "XSw = binarySolubilitySwelling%s(P0,T0,Mp,Ms,**kwargs)" % (theory)
		result = XSw
	
	elif not isListOrNpyArray(T0) and isListOrNpyArray(P0):
		result = [[range(0,len(P0))] for x in range(4)]
		T = range(0,len(P0))
		m_s = range(0,len(P0))
		Sw = range(0,len(P0))
		for i in range(0,len(P0)):
			exec "XSw = binarySolubilitySwelling%s(P0[i],T0,Mp,Ms,**kwargs)" % (theory)
			T[i] = XSw[1]
			m_s[i] = XSw[2]
			Sw[i] = XSw[3]
		result[0] = P0
		result[1] = T
		result[2] = m_s
		result[3] = Sw

	elif not isListOrNpyArray(P0) and isListOrNpyArray(T0):
		result = [[range(0,len(T0))] for x in range(4)]
		P = range(0,len(T0))
		m_s = range(0,len(T0))
		Sw = range(0,len(T0))
		for i in range(0,len(T0)):
			exec "XSw = binarySolubilitySwelling%s(P0,T0[i],Mp,Ms,**kwargs)" % (theory)
			P[i] = XSw[0]
			m_s[i] = XSw[2]
			Sw[i] = XSw[3]
		result[0] = P
		result[1] = T0
		result[2] = m_s
		result[3] = Sw
	
	elif isListOrNpyArray(P0) and isListOrNpyArray(T0):
		result = [[range(0,len(T0))] for x in range(4)]
		P = range(0,len(T0))
		m_s = range(0,len(T0))
		Sw = range(0,len(T0))
		for i in range(0,len(T0)):
			exec "XSw = binarySolubilitySwelling%s(P0[i],T0[i],Mp,Ms,**kwargs)" % (theory)
			P[i] = XSw[0]
			m_s[i] = XSw[2]
			Sw[i] = XSw[3]
		result[0] = P0
		result[1] = T0
		result[2] = m_s
		result[3] = Sw
	
	else:
		raise ValueError('In calculateBinarySwelling: Unknown error involving P0 and T0.')
	
	return result

def calculateBinarySolubility(theory,P0,T0,Mp,Ms,**kwargs):

	if not isListOrNpyArray(P0) and not isListOrNpyArray(T0):
		exec "P,T,m_s,Sw = binarySolubilitySwelling%s(P0,T0,Mp,Ms,**kwargs)" % (theory)
		result = [P,T,m_s]
	
	elif not isListOrNpyArray(T0) and isListOrNpyArray(P0):
		result = [[range(0,len(P0))] for x in range(3)]
		T = range(0,len(P0))
		m_s = range(0,len(P0))
		for i in range(0,len(P0)):
			exec "X = binarySolubilitySwelling%s(P0[i],T0,Mp,Ms,**kwargs)" % (theory)
			T[i] = X[1]
			m_s[i] = X[2]
		result[0] = P0
		result[1] = T
		result[2] = m_s

	elif not isListOrNpyArray(P0) and isListOrNpyArray(T0):
		result = [[range(0,len(T0))] for x in range(3)]
		P = range(0,len(T0))
		m_s = range(0,len(T0))
		for i in range(0,len(T0)):
			exec "X = binarySolubilitySwelling%s(P0,T0[i],Mp,Ms,**kwargs)" % (theory)
			P[i] = X[0]
			m_s[i] = X[2]
		result[0] = P
		result[1] = T0
		result[2] = m_s
	
	elif isListOrNpyArray(P0) and isListOrNpyArray(T0):
		result = [[range(0,len(T0))] for x in range(3)]
		P = range(0,len(T0))
		m_s = range(0,len(T0))
		for i in range(0,len(T0)):
			exec "X = binarySolubilitySwelling%s(P0[i],T0[i],Mp,Ms,**kwargs)" % (theory)
			m_s[i] = X[2]
		result[0] = P0
		result[1] = T0
		result[2] = m_s
	
	else:
		raise ValueError('In calculateBinarySwelling: Unknown error involving P0 and T0.')
	
	return result

def calculateBinarySwelling(theory,P0,T0,Mp,Ms,**kwargs):

	if not isListOrNpyArray(P0) and not isListOrNpyArray(T0):
		exec "P,T,m_s,Sw = binarySolubilitySwelling%s(P0,T0,Mp,Ms,**kwargs)" % (theory)
		result = [P,T,Sw]

	elif not isListOrNpyArray(T0) and isListOrNpyArray(P0):
		result = [[range(0,len(P0))] for x in range(3)]
		T = range(0,len(P0))
		Sw = range(0,len(P0))
		for i in range(0,len(P0)):
			exec "Swlng = binarySolubilitySwelling%s(P0[i],T0,Mp,Ms,**kwargs)" % (theory)
			T[i] = Swlng[1]
			Sw[i] = Swlng[3]
		result[0] = P0
		result[1] = T
		result[2] = Sw
	
	elif not isListOrNpyArray(P0) and isListOrNpyArray(T0):
		result = [[range(0,len(T0))] for x in range(3)]
		P = range(0,len(T0))
		Sw = range(0,len(T0))
		for i in range(0,len(T0)):
			exec "Swlng = binarySolubilitySwelling%s(P0,T0[i],Mp,Ms,**kwargs)" % (theory)
			P[i] = Swlng[0]
			Sw[i] = Swlng[3]
		result[0] = P
		result[1] = T0
		result[2] = Sw
	
	elif isListOrNpyArray(P0) and isListOrNpyArray(T0):
		result = [[range(0,len(T0))] for x in range(3)]
		P = range(0,len(T0))
		Sw = range(0,len(T0))
		for i in range(0,len(T0)):
			exec "Swlng = binarySolubilitySwelling%s(P0[i],T0[i],Mp,Ms,**kwargs)" % (theory)
			Sw[i] = Swlng[3]
		result[0] = P0
		result[1] = T0
		result[2] = Sw
	
	else:
		raise ValueError('In calculateBinarySwelling: Unknown error involving P0 and T0.')
	
	return result

def residualFunction(A0,A,weight=1.0):
	if len(A0) != len(A):
		raise ValueError('In residual: The number of experimental points and number of theoretical points are not equal.')
	
	residual = npy.zeros(len(A0))
	
	for i in range(0,len(A0)):
		residual[i] = weight*((A0[i]-A[i]))#/A0[i]  #Kier original had no hash and no absolute
	
	
	# print 'weight is', weight
	return residual

def binaryResidual(theory,P0_X,T0_X,X0_X,P0_S,T0_S,S0_S,Mp,Ms,fit_type,fs,**kwargs):

	m_s = npy.zeros(len(P0_X))
	Sw = npy.zeros(len(P0_S))
	res_X = npy.zeros(len(P0_X))
	res_S = npy.zeros(len(P0_S))
	
	if len(P0_X) != len(P0_S):
		warnings.warn('In binaryResidual: Mismatch in solubility and swelling data number. Results may be skewed.')

	suppress_print = kwargs.pop('suppress_print',False)
	if not suppress_print:
		for key,value in kwargs.items():
			print '%s=%s' % (key,value)

	if 'X' in fit_type:
		P0,T0,m_s = calculateBinarySolubility(theory,P0_X,T0_X,Mp,Ms,**kwargs)
		print 'solubility weight is'
		res_X = residualFunction(X0_X,m_s,1.0-fs)
	if 'S' in fit_type:
		P0,T0,Sw = calculateBinarySwelling(theory,P0_S,T0_S,Mp,Ms,**kwargs)
		# Sw_dash=npy.array(Sw)-1
		# S0_S_dash=npy.array(S0_S)-1
		print 'swelling weight is'
		res_S = residualFunction(S0_S,Sw,fs)
	
	if 'X' in fit_type and 'S' in fit_type:
		# print 'hurry'
		residual = npy.concatenate((res_X,res_S),axis=0)
	elif 'X' in fit_type:
		residual = res_X
	elif 'S' in fit_type:
		residual = res_S
	else:
		raise ValueError('In binaryResidual: fit_type must contain X and/or S.')

	return residual

def calculateBinaryResidualCHV(params,P0_X,T0_X,X0_X,P0_S,T0_S,S0_S,fit_type,method):

	fs = params['fs'].value
	Mp = params['Mp'].value
	Ms = params['Ms'].value
	vhp = params['vhp'].value
	vhs = params['vhs'].value
	vr = params['vr'].value	
	print vr/vhp

	if 'Ppstar' in params and 'Tpstar' in params and 'Rpstar' in params:
		Ppstar = params['Ppstar'].value
		Tpstar = params['Tpstar'].value
		Rpstar = params['Rpstar'].value
		Tpstar = Tpstar*vr/vhp
		kwargs = {'Ppstar':Ppstar,'Tpstar':Tpstar,'Rpstar':Rpstar}
	elif 'alpha_p' in params and 'vhp' in params and 'epsilon_p' in params:
		alpha_p = params['alpha_p'].value
		vhp = params['vhp'].value
		epsilon_p = params['epsilon_p'].value
		kwargs = {'alpha_p':alpha_p,'vhp':vhp,'epsilon_p':epsilon_p}
	else:
		raise ValueError('In calculateBinaryResidualCHV, pure polymer: Either molecular (alpha_p,vhp,epsilon_p) or characteristic (Ppstar,Tpstar,Rpstar) parameters must be used.')
	
	if 'Psstar' in params and 'Tsstar' in params and 'Rsstar' in params:
		Psstar = params['Psstar'].value
		Tsstar = params['Tsstar'].value
		Rsstar = params['Rsstar'].value
		Tsstar = Tsstar*vr/vhs
		kwargs.update({'Psstar':Psstar,'Tsstar':Tsstar,'Rsstar':Rsstar})
	elif 'alpha_s' in params and 'vhs' in params and 'epsilon_s' in params:
		alpha_s = params['alpha_s'].value
		vhs = params['vhs'].value
		epsilon_s = params['epsilon_s'].value
		kwargs.update({'alpha_s':alpha_s,'vhs':vhs,'epsilon_s':epsilon_s})
	else:
		raise ValueError('In calculateBinaryResidualCHV, pure solvent: Either molecular (alpha_s,vhs,epsilon_s) or characteristic (Psstar,Tsstar,Rsstar) parameters must be used.')
	
	if 'k12' in params and 'delta' in params:
		k12 = params['k12'].value
		delta = params['delta'].value
		kwargs.update({'k12':k12,'delta':delta})
	elif 'zeta' in params and 'delta' in params:
		zeta = params['zeta'].value
		delta = params['delta'].value
		kwargs.update({'zeta':zeta,'delta':delta})
	else:
		raise ValueError('In calculateBinaryResidualCHV, mixture parameters: (k12,delta) or (zeta,delta) mixture parameters must be used.')
	
	if method == 'disparate':
		pass
	elif method == 'single':
		kwargs.update({'method':'single'})
	elif method == 'mixed':
		kwargs.update({'method':'mixed'})
	else:
		warnings.warn("In calculateBinaryResidualCHV, method: method parameter not specified. method = 'disparate' will be used.")
	
	if 'verbose' in params:
		verbose = params['verbose'].value
		kwargs.update({'verbose':verbose})
	
	res = binaryResidual('CHV',P0_X,T0_X,X0_X,P0_S,T0_S,S0_S,Mp,Ms,fit_type,fs,**kwargs)
	
	print '==> Done for values above.'

	return res


if __name__ == "__main__":

	Polymer_Type='PMMA'
	Solvent='CO2'

	Parameters_Paper ='Self_Grassia'			# P*T*R* and g,epsilon_2,x (PVT-Tg Data Paper or Direct P*T*R* Values Reference)
	Cp_Polymer_Weight = '02kilo_POST_THESIS'	# g,epsilon_2,x (Cp Paper Reference)
	Paper_Number = 'Paper15'						# Solubility or Swelling Data Reference
	#for PS: 'Paper4_11_12'
	#for PMMA: 'Paper15'
	kwargs = {'Polymer_Type':Polymer_Type,'Solvent':Solvent,'Parameters_Paper':Parameters_Paper,'Paper_Number':Paper_Number,'Cp_Polymer_Weight':Cp_Polymer_Weight}

	Ppstar,Tpstar,Rpstar,Mp,Psstar,Tsstar,Rsstar,Ms,P_exp,Tg_exp=Parameters_of_Different_Polymers(**kwargs)
	P0_X_complete,T0_X_complete,X0_X_complete,P0_S_complete,T0_S_complete,S0_S_complete,Rubber0_X_complete,Rubber0_S_complete=loadExperimentSwXData(**kwargs)
	Far_Above_Data=False
	P0_X,P0_X_above_Tg,P0_X_far_above_Tg,T0_X,T0_X_above_Tg,T0_X_far_above_Tg,X0_X,X0_X_above_Tg,X0_X_far_above_Tg,Rubber0_X,Rubber0_X_above_Tg,Rubber0_X_far_above_Tg,P0_S,P0_S_above_Tg,P0_S_far_above_Tg,T0_S,T0_S_above_Tg,T0_S_far_above_Tg,S0_S,S0_S_above_Tg,S0_S_far_above_Tg,Rubber0_S,Rubber0_S_above_Tg,Rubber0_S_far_above_Tg = SplitExperimental_X_Sw_Data(P0_X_complete,T0_X_complete,X0_X_complete,P0_S_complete,T0_S_complete,S0_S_complete,Rubber0_X_complete,Rubber0_S_complete,Far_Above_Data,**kwargs)

	number_of_isotherm, result = Split_Isotherms(P0_X,T0_X,X0_X,'X')
	P0_X_T1,T0_X_T1,X0_X_T1,P0_X_T2,T0_X_T2,X0_X_T2,P0_X_T3,T0_X_T3,X0_X_T3,P0_X_T4,T0_X_T4,X0_X_T4,P0_X_T5,T0_X_T5,X0_X_T5,P0_X_T6,T0_X_T6,X0_X_T6,P0_X_T7,T0_X_T7,X0_X_T7,P0_X_T8,T0_X_T8,X0_X_T8 = result[0],result[1],result[2],result[3],result[4],result[5],result[6],result[7],result[8],result[9],result[10],result[11],result[12],result[13],result[14],result[15],result[16],result[17],result[18],result[19],result[20],result[21],result[22],result[23]
	# print P0_X_T1,T0_X_T1,X0_X_T1,P0_X_T2,T0_X_T2,X0_X_T2,P0_X_T3,T0_X_T3,X0_X_T3,P0_X_T4,T0_X_T4,X0_X_T4,P0_X_T5,T0_X_T5,X0_X_T5

	number_of_isotherm_swelling, result = Split_Isotherms(P0_S,T0_S,S0_S,'S')
	P0_S_T1,T0_S_T1,S0_S_T1,P0_S_T2,T0_S_T2,S0_S_T2,P0_S_T3,T0_S_T3,S0_S_T3,P0_S_T4,T0_S_T4,S0_S_T4,P0_S_T5,T0_S_T5,S0_S_T5,P0_S_T6,T0_S_T6,S0_S_T6,P0_S_T7,T0_S_T7,S0_S_T7,P0_S_T8,T0_S_T8,S0_S_T8 = result[0],result[1],result[2],result[3],result[4],result[5],result[6],result[7],result[8],result[9],result[10],result[11],result[12],result[13],result[14],result[15],result[16],result[17],result[18],result[19],result[20],result[21],result[22],result[23]
	# print P0_S_T1,T0_S_T1,S0_S_T1,P0_S_T2,T0_S_T2,S0_S_T2,P0_S_T3,T0_S_T3,S0_S_T3,P0_S_T4,T0_S_T4,S0_S_T4,P0_S_T5,T0_S_T5,S0_S_T5

	vhp = kB*Tpstar/Ppstar
	vhs = kB*Tsstar/Psstar
	vr = vhp
	# vr = (vhp+vhs)/2
	# vr = math.sqrt(vhp*vhs)
	# Tsstar = Tsstar*vr/vhs
	# Tpstar = Tpstar*vr/vhp

	# P0_X = P0_X_T2
	# T0_X = T0_X_T2
	# X0_X = X0_X_T2

	# P0_S = P0_S_T2
	# T0_S = T0_S_T2
	# S0_S = S0_S_T2

	print 'P0_X=',P0_X
	print 'T0_X=', T0_X
	print 'X0_X=', X0_X
	print 'Rubber0_X=', Rubber0_X

	print 'P0_S=',P0_S
	print 'T0_S=', T0_S
	print 'S0_S=', S0_S
	print 'Rubber0_S=', Rubber0_S

	# P0_X = npy.concatenate((P0_X,[0.0001]),axis=0)
	# T0_X = npy.concatenate((T0_X,[373.15]),axis=0)
	# X0_X = npy.concatenate((X0_X,[0]),axis=0)

	#Initializing the parameters.
	params = Parameters()
	#The following code sets up the model's parameters. It includes both fitting parameters and parameters that will remain fixed
	#for the fitting. The values given are the inital guesses of fitting parameters and values of fixed parameters.
	#           	(Name,		Value,		Vary?,	Min,	Max,	Expr)
	params.add_many(('zeta',	1.03,		True,	None,	None,	None),
					('delta',	1.41,		True,	None,	None,	None),
					('Ppstar',	Ppstar,		False,	None,	None,	None),
					('Tpstar',	Tpstar,		False,	None,	None,	None),
					('Rpstar',	Rpstar,		False,	None,	None,	None),
					('Mp',		Mp,			False,	None,	None,	None),
					('Psstar',	Psstar,		False,	None,	None,	None),
					('Tsstar',	Tsstar,		False,	None,	None,	None),
					('Rsstar',	Rsstar,		False,	None,	None,	None),
					('Ms',		Ms,			False,	None,	None,	None),
					('vhs',		vhs,		False,	None,	None,	None),
					('vhp',		vhp,		False,	None,	None,	None),
					('vr',		vr,			False,	None,	None,	None),
					('fs',		0.5,		False,	None,	None,	None),
					('verbose',	False,		False,	None,	None,	None))
	#fs=0.0 => Pure X0 fit.
	#fs=Weight of swelling residual in simultaneous solubility-swelling fit
	#Running the Levenberg-Marquart algorithm on the residuals in order to do least squares fitting. This will return the fitted value of the RESIDUALS.
	#These need to be added to the experimental datapints to find the fitted pressures.
	print('For POLYMER: {} and SOLVENT: {}.'.format(Polymer_Type,Solvent))
	print('Using {} parameters.'.format(Parameters_Paper))

	fit = minimize(calculateBinaryResidualCHV,params,args=(P0_X,T0_X,X0_X,P0_S,T0_S,S0_S,['X','S'],'disparate'))

	#Reporting the values of the parameters. NEED TO FIGURE OUT HOW TO PRINT THIS TO FILE.
	report_fit(fit.params)
