import numpy as np 
import HW2_Functions

find_deltatstar = HW2_Functions.find_deltatstar

'''
Question 1:
Determine the largest value of h, for which the AB3 method is absolutely stable
'''
A = np.array([[0, 10, -10],
              [-100, -1, 0],
              [0, 10, -100]])

def AB3_boundary(theta):
  y = (12*np.e**(3j*theta) - 12*np.e**(2j*theta))/(5-16*np.e**(theta*1j)+23*np.e**(2j*theta)) 
  return y
delta_t_star = find_deltatstar(A, AB3_boundary, "Three-step Adams-Bashforth", 1)


