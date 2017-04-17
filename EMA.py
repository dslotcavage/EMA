import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

#Source for Maxwell Garnett Forumla: wikipedia
def MaxwellGarnettFormula(eps_base,eps_incl,vol_incl):
	small_number_cutoff = 1e-6

	if vol_incl < 0 or vol_incl > 1:
		print('WARNING: volume portion of inclusion material is out of range!')
	factor_up = 2*(1-vol_incl)*eps_base+(1+2*vol_incl)*eps_incl
	factor_down = (2+vol_incl)*eps_base+(1-vol_incl)*eps_incl
	if min(abs(factor_down)) < small_number_cutoff:
		print('WARNING: the effective medium is singular!')
		eps_mean = 0
	else:
		eps_mean = np.multiply(eps_base,np.divide(factor_up,factor_down))

	return(eps_mean)

#Source https://www.physics.ohio-state.edu/~stroud/stroud1.pdf Eq. 1
#The effective medium approximations: Some recent developments
#David Stroud
#Variables use convention from the source, ea is epislon for material a, p is volume fraction
def BruggemanForumula(ea,eb,p):
	operator = 8*np.multiply(ea,eb) + np.square(ea - 2*eb - 3*ea*p + 3*eb*p)
	eps_effective = 0.25*(-ea + 2*eb +3*ea*p - 3*eb*p - np.sqrt(operator))
	return(eps_effective)

#Hummel 10.13 and 10.14
def nk_to_permittivity(n,k):
	eps = np.multiply(n,n) - np.multiply(k,k)
	eps = eps-1j*2 * np.multiply(n,k)
	return(eps)

#Hummel 10.15 and 10.16
def permittivity_to_nk(eps):
	n = np.sqrt(0.5*(np.sqrt(np.square(eps.real)+np.square(eps.imag))+eps.real))
	k = np.sqrt(0.5*(np.sqrt(np.square(eps.real)+np.square(eps.imag))-eps.real))
	return(n,k)

def import_nk(filename):
	nk = np.genfromtxt(filename,delimiter = ',', skip_header = 1)
	return(nk)

def nk_smoothing(m_nk,f_nk):

	#determining the inner bounds for lambda (so as not to extrapolate raw data)
	if m_nk[0,0] < f_nk[0,0]:
		lambda_min = f_nk[0,0]
	else:
		lambda_min = m_nk[0,0]

	if m_nk[len(m_nk)-1,0] > f_nk[len(f_nk)-1,0]:
		lambda_max = f_nk[len(f_nk)-1,0]
	else:
		lambda_max = m_nk[len(m_nk)-1,0]

	lambdanew = np.linspace(lambda_min, lambda_max, num = ((lambda_max-lambda_min)*1000+1))
 	m_n_interp = interp1d(m_nk[:,0],m_nk[:,1], kind = 'cubic')
 	m_k_interp = interp1d(m_nk[:,0],m_nk[:,2], kind = 'cubic')
 	f_n_interp = interp1d(f_nk[:,0],f_nk[:,1], kind = 'cubic')
 	f_k_interp = interp1d(f_nk[:,0],f_nk[:,2], kind = 'cubic')

	return(m_n_interp, m_k_interp, f_n_interp, f_k_interp, lambdanew)

# Name of the nk csv file (taken from the formate from RefractiveIndex.info)
metal_file = 'Cu_nk.csv'
fluid_file = 'H2O_nk.csv'

# Volume fraction of metal inclusion
f = 0.1

#Column 0 = wavelength (um), col 1 = n, col 2 = k
metal_nk = import_nk(metal_file)
fluid_nk = import_nk(fluid_file)

#interpolates data and returns functions such that data of different lengths can be used for MaxwellGarnett
metal_n, metal_k, fluid_n, fluid_k, wavelength = nk_smoothing(metal_nk,fluid_nk)

metal_eps = nk_to_permittivity(metal_n(wavelength),metal_k(wavelength))
fluid_eps = nk_to_permittivity(fluid_n(wavelength),fluid_k(wavelength))



# caluclations based on Bruggeman model
medium_eps = BruggemanForumula(fluid_eps,metal_eps,f)
medium_n, medium_k = permittivity_to_nk(medium_eps)

for i in range(6,10):
	f=0.1*i
	medium_eps = BruggemanForumula(fluid_eps,metal_eps,f)
	plt.plot(wavelength,medium_eps.imag,'.')
plt.show()

# refl_numerator = np.square(medium_n-1)+np.square(medium_k)
# refl_denominator = np.square(medium_n+1)+np.square(medium_k)
# refl = np.divide(refl_numerator,refl_denominator)

# #import data for comparison
# ITO_refl = np.genfromtxt('10Ohmsq_reflection.txt',skip_header = 19, skip_footer = 35)
# PbCu_refl = np.genfromtxt('PbCureflection.txt',skip_header = 19, skip_footer = 35)
# #converting from nm to um and % to fraction
# ITO_refl[:,0] = .001*ITO_refl[:,0]
# PbCu_refl[:,0] = .001*PbCu_refl[:,0]
# ITO_refl[:,1] = .01*ITO_refl[:,1]
# PbCu_refl[:,1] = .01*PbCu_refl[:,1]
# #interpolating to be able to use the same wavelength range
# sim_refl = interp1d(wavelength,refl)

# plt.plot(wavelength,refl)
# plt.plot(PbCu_refl[:,0],PbCu_refl[:,1])
# plt.ylim(0,1)	
# plt.show()

# #saving as a CSV to use in transfer matrix modeling
# EMAoutput = np.column_stack((1000*wavelength, medium_n, medium_k))
# np.savetxt('nk_EMA.csv', EMAoutput, delimiter = ',', header = 'Wavelength, n, k')

# # calculations based on numerical weighted average of the two materials' optical constants
# for i in range(0,10):
# 	f = 0.1*i
# 	effective_n = f*metal_n(wavelength)+(1-f)*fluid_n(wavelength)
# 	effective_k = f*metal_k(wavelength)+(1-f)*fluid_k(wavelength)
# 	refl_numerator = np.square(effective_n-1)+np.square(effective_k)
#  	refl_denominator = np.square(effective_n+1)+np.square(effective_k)
# 	refl = np.divide(refl_numerator,refl_denominator)

# 	plt.plot(wavelength,refl)
# plt.show()



# # caluclations based on MaxwellGarnett model
# for i in range (0,10):
# 	f = 0.1*i
# 	medium_eps1 = MaxwellGarnettFormula(fluid_eps1,metal_eps1,f)
# 	medium_eps2 = MaxwellGarnettFormula(fluid_eps2,metal_eps2,f)

# 	medium_n, medium_k = permittivity_to_nk(medium_eps1,medium_eps2)

# 	refl_numerator = np.square(medium_n-1)+np.square(medium_k)
# 	refl_denominator = np.square(medium_n+1)+np.square(medium_k)
# 	refl = np.divide(refl_numerator,refl_denominator)

# 	plt.plot(wavelength,refl)
# plt.show()

