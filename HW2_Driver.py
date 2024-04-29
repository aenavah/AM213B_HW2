import numpy as np 
import matplotlib.pyplot as plt
import HW2_Functions

find_deltatstar = HW2_Functions.find_deltatstar
plot_boundary = HW2_Functions.plot_boundary
AB3 = HW2_Functions.AB3
RK3 = HW2_Functions.RK3
RK3_absolute_stability_region = HW2_Functions.RK3_absolute_stability_region

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
delta_t_star = find_deltatstar(A, AB3_boundary, "Three-step Adams-Bashforth")


'''
Question 2b
plot the boundary of absolute stability region for BDF3
'''
plt.clf()
def BDF3_boundary(theta):
  y = (11*np.e**(3j*theta) - 18*np.e**(2j*theta) + 9*np.e**(1j*theta) - 2)/(6*np.e**(3j*theta))
  return y
plot_boundary(BDF3_boundary)
plt.title("BDF3 Absolute Stability Region")
plt.savefig("BDF_Absolute_Stability_Region.jpg")

'''
Question 1c
'''
for dt in [10**(-4), delta_t_star + 0.0001, delta_t_star - 0.0001]:
  T = 10
  n = int(T/dt)
  ts = np.arange(0, 10, dt)
  ts = ts[0 : n]
  u = np.zeros((n, 3))

  u[0] = np.array([10, 10, 10]) 

  #call RK method from Hw1
  for i in range(1, 2+1):
    u_next = RK3(A, u[i-1], dt)
    u[i] = u_next
    u_last = u_next

  #calling next iterations from AB3 method
  for i in range(3, n):
    u_2 = u[i - 1]
    u_1 = u[i - 2]
    u_0 = u[i - 3]

    u_3 = AB3(A, u_0, u_1, u_2, dt)
    u[i] = u_3
  
  dy1 = u[:, 0]
  dy2 = u[:, 1]
  dy3 = u[:, 2]

  plt.clf()
  plt.plot(ts, dy1, label = "dy_1")
  plt.plot(ts, dy2, label = "dy_2")
  plt.plot(ts, dy3, label = "dy_3")
  plt.legend()
  plt.title("delta_t = " + str(dt))
  plt.savefig("AM3_Solution_dt=" + str(dt) + ".jpg")

plt.clf()
'''Q4b'''
A = np.array([[0, 0, 0], 
              [1/2, 0, 0], 
              [0, 3/4, 0]])
b = np.array([[2/9],[1/3],[4/9]])
h = np.ones_like(b)
print(b.shape)
print(h.shape)
RK3_absolute_stability_region(A, b, h)

  