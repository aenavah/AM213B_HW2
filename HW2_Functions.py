import numpy as np
import matplotlib.pyplot as plt
from math import atan2


def BDF3_boundary(theta):
  y = (12*np.e**(3j*theta) - 12*np.e**(2j*theta))/(5-16*np.e**(theta*1j)+23*np.e**(2j*theta)) 
  return y

def q1a(matrix):
  eigenvalues = np.linalg.eigvals(matrix)
  thetas = np.linspace(0, 2*np.pi)

  real_vals = []
  imag_vals = []
  
  #plotting actual boundary
  for theta in thetas:
    y = BDF3_boundary(theta) 
    real_part = y.real
    imag_part = y.imag
    real_vals.append(real_part)
    imag_vals.append(imag_part)
  
  plt.grid()
  plt.xlabel("Re(z)")
  plt.ylabel("Im(z)")
  plt.plot(real_vals, imag_vals)  
  for delta_t in np.arange(0, 1, 10**(-5)):
    for lambda_j in eigenvalues:
    
      #plotting actual boundary
      lambda_scaled = lambda_j * delta_t

      lambda_real = lambda_scaled.real
      lambda_imag = lambda_scaled.imag 
      
      #getting angle of lambda for comparison:
      #angle = atan2(lambda_imag, lambda_real)
      #comp = BDF3_boundary(angle)
      
      stable = 0
      #iterate through all x and y values and check if lambda less than all for each delta t 
      for i in range(len(real_vals)):
        if (abs(lambda_real) < real_vals[i]) and (abs(lambda_imag) < imag_vals[i]):
          stable = 1
          plt.plot([lambda_real], [lambda_imag], 'ro') 
          break
  plt.show()
      

  #for lambda_j in eigenvalues:
    #for delta_t in range(0,1,10**(-3)):
    #delta_t = 10**(-3)

  
