import numpy as np 
import HW2_Functions

q1a = HW2_Functions.q1a

'''
Question 1:
Determine the largest value of h, for which the AB3 method is absolutely stable
'''
A = np.array([[0, 10, -10],
              [-100, -1, 0],
              [0, 10, -100]])

q1a(A)


