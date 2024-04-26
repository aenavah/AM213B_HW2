import numpy as np
import matplotlib.pyplot as plt
from math import atan2


def find_deltatstar(matrix, method, method_name):
  eigenvalues = np.linalg.eigvals(matrix)
  thetas = np.linspace(0, 2*np.pi)

  real_vals = []
  imag_vals = []
  
  #plotting actual boundary
  for theta in thetas:
    y = method(theta) 
    real_part = y.real
    imag_part = y.imag
    real_vals.append(real_part)
    imag_vals.append(imag_part)
  plt.grid()
  plt.xlabel("Re(z)")
  plt.ylabel("Im(z)")
  plt.plot(real_vals, imag_vals)  

  #iterating through delta_ts 
  for delta_t in np.arange(1, 0, -10**(-5)):

    lambda_scaled_reals = []
    lambda_scaled_imags = []

    for lambda_j in eigenvalues:

      #getting lambda at this scale 
      lambda_scaled = lambda_j * delta_t
      lambda_real = lambda_scaled.real
      lambda_imag = lambda_scaled.imag 
      
      #getting angle of lambda for comparison:
      angle = atan2(lambda_imag, lambda_real)
      bound_at_angle = method(angle)
      realbound_at_angle = bound_at_angle.real
      imagbound_at_angle = bound_at_angle.imag

      #if lambda is found for this delta_t:
      if (abs(realbound_at_angle) > abs(lambda_real)) and (abs(imagbound_at_angle) > abs(lambda_imag)):
        lambda_scaled_reals.append(lambda_real)
        lambda_scaled_imags.append(lambda_imag)
        continue

    #if all eigenvalues are in region:
    if len(lambda_scaled_reals) == len(eigenvalues):
      print(lambda_scaled_reals)
      print(delta_t)
      plt.scatter(lambda_scaled_reals, lambda_scaled_imags, marker = "o", color = "red", facecolors = "none")
      plt.title(method_name + " with delta_t = " + str(delta_t) )
      plt.show()
      return




      


  
