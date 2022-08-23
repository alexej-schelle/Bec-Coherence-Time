############################################################################################################################################################################
#
#
# -*- coding: utf-8 -*-
#
#
############################################################################################################################################################################
#                                                                                                                                              
#																										  
#
#   Markov sampling method (generalised Metropolis-Hastings sampling algorithm) for calculating the
#   coherence time (distribution) of the order parameter for a given total particle number, 
#   temperature and trap frequency of a Bose-Einstein condensate.
#
#   This software version corresponds to publication of Fluctuation and Noise Letters, Vol. 16, No. 01, 1750009 (2017).
#
#  * :
# 
#   License Copyright:  Dr. A. Schelle, Sudetenstr. 76, 87600 Kaufbeuren 
#   License Type :      MIT license (2017)
#   License Contact:    E-Mail : alexej.schelle@gmail.com
#
#   Free support is provided per email at support@krealix.de.   
#   For more explicit consulting and discussions, please schedule a meeting at https://calendly.com/alexej-schelle/.  
# 
#   ** : 
#
#   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files 
#   (the "Software" bec_symmetry_breaking.py), to deal in the Software without restriction, including without limitation the rights to use, 
#   copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
#   furnished to do so, subject to the following conditions:
# 
#   The above copyright notice (*) and this permission notice (**) shall be included in all copies or substantial portions of the Software.
# 
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
#   WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#   In order to become sponsor of the project, simply get in touch with krealix.de at support@krealix.de.
#
############################################################################################################################################################################

import os
import sys
import math
import random
import numpy
import numpy as np
import pylab
import matplotlib.pyplot as plt
import operator

maxmode = 500 # Typical mode size for analysis : 50 - 2500 modes.
ptn = 1000 # Typical particle number : 10^3 - 10^5
sample = 1000 # Typical sample size : 10^5 - 10^8

omx = 2.0*math.pi*42.0 # Trap frequency in x direction
omy = 2.0*math.pi*42.0 # Trap frequency in y direction 
omz = 2.0*math.pi*120.0 # Trap frequency in z direction

# omx = 2.0*math.pi*50.0 # Trap frequency in x direction
# omy = 2.0*math.pi*75.0 # Trap frequency in y direction 
# omz = 2.0*math.pi*125.0 # Trap frequency in z direction

start_temp = 10.0 # in units of nK

rhbkb = 7.63822291E-3 # in units of nK
norm = 0.0
int_norm = 0
drop = 0
prob = 0.0
phase_0 = 0.0
z_start = 0.0
mu_start = 0.0

print 'Trap depth [nK]: ', maxmode*omz*rhbkb # around few muK
print 'Critical temperature [nK]: ',  rhbkb*pow(ptn,1.0/3.0)*pow(omx*omy*omz,1.0/3.0)/pow(1.202,1.0/3.0) # around 1.15 mK

en_x = ['']*maxmode
en_y = ['']*maxmode
en_z = ['']*maxmode

pols_x = ['']*maxmode
pols_y = ['']*maxmode
pols_z = ['']*maxmode

pols = ['']*maxmode

x = ['']*ptn
y = ['']*ptn

mu_x = [] # Collect the real part of two-dimensional complex time.
mu_y = [] # Collect the imaginary part of two-dimensional complex time.
mu_abs = [] # Collect the absolute value of two-dimensional complex time.

for l in range(1, sample):
          
    drop = 1
    z = random.uniform(1, ptn) # Random sample of non-condensate particle number
    temp = start_temp
    
    mu = 0.0
    norm = 0.0
        
    print 'Sample step Nr. ' + str(l)
                  
    for k in range(1, maxmode):
        
        en_x[k] = rhbkb*k*omx/temp # Energy in x direction
        en_y[k] = rhbkb*k*omy/temp # Energy in y direction
        en_z[k] = rhbkb*k*omz/temp # Energy in z direction
     
        pols_x[k] = 1.0/(1.0-math.exp(-en_x[k])) # Complex poles in x direction
        pols_y[k] = 1.0/(1.0-math.exp(-en_y[k])) # Complex poles in y direction
        pols_z[k] = 1.0/(1.0-math.exp(-en_z[k])) # Complex poles in z direction
                    
    for k in range(1, maxmode):
    
        pols[k-1] = pols_x[maxmode-k]*pols_y[maxmode-k]*pols_z[maxmode-k]-1.0 # General poles
    
    pols[maxmode-1] = -(ptn-z_start)
    
    prob = complex(0.0,0.0)
    p = ['']*maxmode
    phase_0 = 0.0
   
    x = numpy.roots(pols) # Complex roots of the number conserving equation
    
    for k in range(0,len(x)):

        p[k] = random.uniform(0.0,1.0)   # Random amplitudes - set p[k] = delta(k-k_random) for single BEC spectrum
        norm = norm + p[k]*p[k]*(x[k].real**2+x[k].imag**2) # Total norm
    
    norm = math.sqrt(norm)
    
    for k in range(0,len(x)): # Calculate phase of the condensate wave field
        	
        p[k] = p[k]/norm # Random amplitudes
    	
    	if (operator.gt(x[k].real**2 + x[k].imag**2,0.0)):
                
            if (operator.gt(x[k].real,0.0)):
            
                phase_0 = math.atan(x[k].imag/x[k].real)
        
    	    if (operator.iand(operator.lt(x[k].real,0.0),operator.ge(x[k].imag,0.0))):
            
                phase_0 = math.atan(x[k].imag/x[k].real) + math.pi
    
            if (operator.iand(operator.lt(x[k].real,0.0),operator.lt(x[k].imag,0.0))):
	            
                phase_0 = math.atan(x[k].imag/x[k].real) - math.pi
    
            if (operator.iand(operator.eq(x[k].real,0.0),operator.gt(x[k].imag,0.0))):
            
                phase_0 = 0.5*math.pi

            if (operator.iand(operator.eq(x[k].real,0.0),operator.lt(x[k].imag,0.0))):
            
                phase_0 = -0.5*math.pi
        	       
    	mu += x[k]*p[k] # Random amplitudes times phases
        
        if operator.ne(phase_0,0.0):
    	
    	    prob += complex(0.5*p[k]*p[k]*math.log(math.fabs(x[k].real**2 + x[k].imag**2)), p[k]*p[k]*phase_0) # Calculate transition probability
    	
        mu_x.append(prob.real/(omx+omy+omz)/0.5*1.0E3)
        mu_y.append(prob.imag/(omx+omy+omz)/0.5*1.0E3)

        prob = math.sqrt(prob.real**2 + prob.imag**2)

    if (operator.gt(min((math.exp(prob))/(math.exp(mu_start)),1.00),random.uniform(0.00,1.00))): # Condition for transition to another state at equilibrium
 
        mu_start = prob
        z_start = z
        drop = 0
        
        mu_abs.append(prob/(omx+omy+omz)/0.5*1.0E3)
        
plt.figure(1)
plt.hist(mu_x, bins = 500, normed = True)
plt.tick_params(axis='both', which='major', labelsize = 14)
plt.xlabel('$t_1 [ms] $', fontsize = 18)
plt.ylabel('$\pi[t_1]$', fontsize = 18)
plt.savefig('/Users/dr.a.schelle/Desktop/coherence_time_publication/pub_bec_coherence_time/bec_time_distribution/fig_field_modes_1.png')

plt.figure(2)
plt.hist(mu_y, bins = 500, normed = True)
plt.tick_params(axis='both', which='major', labelsize = 14)
plt.xlabel('$t_2 [ms] $', fontsize = 18)
plt.ylabel('$\pi[t_2]$', fontsize = 18)
plt.savefig('/Users/dr.a.schelle/Desktop/coherence_time_publication/pub_bec_coherence_time/bec_time_distribution/fig_field_modes_2.png')

plt.figure(3)
plt.hist(mu_abs, bins = 500, normed = True)
plt.tick_params(axis='both', which='major', labelsize = 14)
plt.xlabel('$t [ms] $', fontsize = 18)
plt.ylabel('$\pi[t]$', fontsize = 18)
plt.savefig('/Users/dr.a.schelle/Desktop/coherence_time_publication/pub_bec_coherence_time/bec_time_distribution/fig_field_modes_3.png')
