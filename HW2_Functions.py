import numpy as np
import matplotlib.pyplot as plt
from math import atan2


def plot_boundary(method):
  #plots boundary given the boundary function
  # note: need to plt.title and plt.show after calling this function
  plt.clf()
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


def find_deltatstar(matrix, method, method_name):
#finds delta t star given the ODE and the absolute stability boundary
  eigenvalues = np.linalg.eigvals(matrix)
  plot_boundary(method)

  #iterating through delta_ts 
  tolerance = -10**(-4)
  for delta_t in np.arange(1, 0, tolerance):

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
      #print(lambda_scaled_reals)
      #print(delta_t)
      plt.scatter(lambda_scaled_reals, lambda_scaled_imags, marker = "o", color = "red", facecolors = "none", label = "Eigenvalues")
      plt.legend()
      plt.title(method_name + " with delta_t = " + str(delta_t) )
      plt.savefig("scaled_eigenvalues_tol=" + str(tolerance) + ".jpg")
      return delta_t

def f(A, y):
  f = A @ y
  return f 

def RK3(A, y0, dt): 
    '''RK method from HW1'''
    K1 = f(A, y0)
    K2 = f(A, y0 + ((1.0/3.0) * dt * K1))
    K3 = f(A, y0 + ((2.0/3.0) * dt * K2))
    y_next = y0 + (dt * (1.0/4.0)) * (K1 + 3.0 * K3)
    return y_next

def AB3(A, u_k, u_k1, u_k2, delta_t):    
  u_k3 = u_k2 + delta_t*(1/12)*(23*f(A, u_k2) - 16*f(A, u_k1) + 5*f(A, u_k))
  return u_k3
  

def RK3_S(A, b, h, x, y):
  z = complex(x,y)
  S = abs(np.linalg.det(np.eye(3) - (z*A) + z*(h@(b.T))))
  return S

def RK3_absolute_stability_region(A, b, h):
  nx, ny = (200, 200)
  imag = np.linspace(-4, 4, nx)
  real = np.linspace(-4, 4, ny)
  real_v, imag_v = np.meshgrid(real, imag)

  S_grid = np.ones((nx, ny))

  for x_i in range(nx):
    for y_i in range(ny):
      boundary = RK3_S(A, b, h, real_v[x_i, y_i], imag_v[x_i, y_i])
      S_grid[x_i, y_i] = boundary 

  plt.xlabel("Re(z)")
  plt.ylabel("Im(z)")
  plt.contour(real_v, imag_v, S_grid, [1])
  plt.grid()
  plt.title("Boundary of RK3 Method")
  plt.savefig("RK3_Boundary.jpg")
