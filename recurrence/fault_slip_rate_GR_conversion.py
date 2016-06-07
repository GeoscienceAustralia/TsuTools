"""Script to convert between fault slip-rates and Gutenberg-Richter a and b values
Jonathan Griffin, AIFDR
October 2012
"""

import math
import  numpy as np

# Recurrence parameters calculated from fitting seimicity data to a G-R curve
M_min = 4.9 # Minimum magnitude (for catalogue copmmleteness and G-R generation
M_max = 7.0 # Maximum magnitude in source zone
lambda_min =  1.34 #0.175 ##2.925 Number of earthquakes greater than M_min per year (= A_min)
b = 0.94 # Gutenberg-Richter b-value
beta = 2.303*b

# Fault Parameter
slip_fault = 30.0 # Reported slip rate for the fault (mm/a)
L = 370 #56. #146 #93 #299    # fault length (km)
W = 45. # 18. # 14.8 #18. .#18 # fault width
A = L*W # Fault area
A = A*np.power(10,10) # convert to cm^2


slip_fault = slip_fault/10 # convert to cm/a

mu = 3e11 # shear modulus dyn/cm^2
#mu = 3e16 # shear modulus N/km^2

c = 1.5 # Hanks and Kanamori coefficient

# Convert moment magnitude to moment for M_max
moment_max = np.power(10,(c*M_max + 16.1)) # dyn/cm = N/km

slip_seismic = (b*lambda_min*moment_max*np.exp(-1*beta*(M_max-M_min)))/(mu*A*(c-b)*(1-np.exp(-1*beta*(M_max-M_min))))
moment_seismic = (b*lambda_min*moment_max*np.exp(-1*beta*(M_max-M_min)))/((c-b)*(1-np.exp(-1*beta*(M_max-M_min))))

print 'Slip rate derived from seismic activity rate (mm)', slip_seismic*10
print 'Slip rate from reported fault data (mm)', slip_fault*10

# Calculate activity rate from slip rate as check
lambda_min_slip = (mu*A*(slip_fault)*(c-b)*(1-np.exp(-1*beta*(M_max-M_min))))/(b*moment_max*np.exp(-1*beta*(M_max-M_min)))

# Momement rate from slip rate
moment_slip = mu*A*slip_fault

print 'Activity rate (lambda_min) derived fault slip rate', lambda_min_slip
print 'Activity rate (lambda_min) from seismicity data', lambda_min

print 'Moment rate from slip-rate', moment_slip
print 'Moment rate from seismicity', moment_seismic

