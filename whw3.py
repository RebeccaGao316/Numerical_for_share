# This is a sample Python script.
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special
import math
from scipy.integrate import quad



# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

# x = np.linspace(-1,1,100)
# fig,ax = plt.subplots()
# y0 = math.sqrt(2.0)/2 + 0*x
# ax.plot(x,y0,label = '0')
#
# y1 = math.sqrt(6.0)/2 *x
# ax.plot(x,y1,label = '1')
#
# y2 = math.sqrt(45.0/8) *(x**2-1/3)
# ax.plot(x,y2,label = '2')
#
# y3 = math.sqrt(175.0/8) *(x**3-3*x/5)
# ax.plot(x,y3,label = '3')
#
# y4 = math.sqrt(11025.0/128) *(x**4-6*x**2/7+3.0/35)
# ax.plot(x,y4,label = '4')
#
# y5 = math.sqrt(43659.0/128) *(x**5-9*x**3/10+11*x/63)
# ax.plot(x,y5,label = '5')
#
# plt.xlim(-1,1)
# plt.ylim(-1.2,1.2)
#
# plt.show()
def f_cos(x):
    #return np.cos(x)
    return np.sin(x)
def y_0(x):
    return math.sqrt(2.0)/2
def y_1(x):
    return math.sqrt(6.0)/2 *x
def y_2(x):
    return math.sqrt(45.0/8) *(x**2-1/3)
def y_3(x):
    return  math.sqrt(175.0/8) *(x**3-3*x/5)
def y_4(x):
    return math.sqrt(11025.0/128) *(x**4-6*x**2/7+3.0/35)
def y_5(x):
    return math.sqrt(43659.0/128) *(x**5-10*x**3/9+15*x/63)

#first function: f = cos(x)
c = [0,0,0,0,0,0]
c[0] = quad(lambda x: f_cos(x)*y_0(x),-1.0,1.0)[0]

c[1] = quad(lambda x: f_cos(x)*y_1(x),-1.0,1.0)[0]

c[2] = quad(lambda x: f_cos(x)*y_2(x),-1.0,1.0)[0]

c[3] = quad(lambda x: f_cos(x)*y_3(x),-1.0,1.0)[0]


c[4] = quad(lambda x: f_cos(x)*y_4(x),-1.0,1.0)[0]


c[5] = quad(lambda x: f_cos(x)*y_5(x),-1.0,1.0)[0]
print(c[5])

x = np.linspace(-1,1,100)
fig,ax = plt.subplots()
#cos fig1
y_cos_approxi = c[0]*y_0(x)+c[1]*y_1(x)+c[2]*y_2(x)+c[3]*y_3(x)+c[4]*y_4(x)+c[5]*y_5(x)
ax.plot(x,y_cos_approxi,label = 'y_sin_approxi_k=5')
y_cos_approxi = c[0]*y_0(x)+c[1]*y_1(x)+c[2]*y_2(x)+c[3]*y_3(x)+c[4]*y_4(x)
ax.plot(x,y_cos_approxi,label = 'y_sin_approxi_k=4')
y_cos_approxi = c[0]*y_0(x)+c[1]*y_1(x)+c[2]*y_2(x)+c[3]*y_3(x)
ax.plot(x,y_cos_approxi,label = 'y_sin_approxi_k=3')
y_cos_approxi = c[0]*y_0(x)+c[1]*y_1(x)+c[2]*y_2(x)
ax.plot(x,y_cos_approxi,label = 'y_sin_approxi_k=2')
y_cos_approxi = c[0]*y_0(x)+c[1]*y_1(x)
ax.plot(x,y_cos_approxi,label = 'y_sin_approxi_k=1')
y_cos_approxi = c[0]*y_0(x)+0*x
ax.plot(x,y_cos_approxi,label = 'y_sin_approxi_k=0')
y_cos_real = f_cos(x)
ax.plot(x,y_cos_real,label = 'y_sin_real')
leg = plt.legend(loc='lower center')
plt.show()
# # #cos fig2
# y_cos_approxi_dif = c[0]*y_0(x)+c[1]*y_1(x)+c[2]*y_2(x)+c[3]*y_3(x)+c[4]*y_4(x)+c[5]*y_5(x)-f_cos(x)
# ax.plot(x,y_cos_approxi_dif,label = 'y_sin_approxi_dif=5')
# y_cos_approxi_dif = c[0]*y_0(x)+c[1]*y_1(x)+c[2]*y_2(x)+c[3]*y_3(x)+c[4]*y_4(x)-f_cos(x)
# ax.plot(x,y_cos_approxi_dif,label = 'y_sin_approxi_dif=4')
# y_cos_approxi_dif = c[0]*y_0(x)+c[1]*y_1(x)+c[2]*y_2(x)+c[3]*y_3(x)-f_cos(x)
# ax.plot(x,y_cos_approxi_dif,label = 'y_sin_approxi_dif=3')
# y_cos_approxi_dif = c[0]*y_0(x)+c[1]*y_1(x)+c[2]*y_2(x)-f_cos(x)
# ax.plot(x,y_cos_approxi_dif,label = 'y_sin_approxi_dif=2')
# y_cos_approxi_dif = c[0]*y_0(x)+c[1]*y_1(x)-f_cos(x)
# ax.plot(x,y_cos_approxi_dif,label = 'y_sin_approxi_dif=1')
# y_cos_approxi_dif = c[0]*y_0(x)+0*x-f_cos(x)
# ax.plot(x,y_cos_approxi_dif,label = 'y_sin_approxi_dif=0')
# leg = plt.legend(loc='upper center')
# plt.show()
#cos p3
# error_5 = quad(lambda x:(c[0]*y_0(x)+c[1]*y_1(x)+c[2]*y_2(x)+c[3]*y_3(x)+c[4]*y_4(x)+c[5]*y_5(x)-f_cos(x))**2,-1.0,1.0)[0]/quad(lambda x:f_cos(x)**2,-1.0,1.0)[0]
# error_4 = quad(lambda x:(c[0]*y_0(x)+c[1]*y_1(x)+c[2]*y_2(x)+c[3]*y_3(x)+c[4]*y_4(x)-f_cos(x))**2,-1.0,1.0)[0]/quad(lambda x:f_cos(x)**2,-1.0,1.0)[0]
# error_3 = quad(lambda x:(c[0]*y_0(x)+c[1]*y_1(x)+c[2]*y_2(x)+c[3]*y_3(x)-f_cos(x))**2,-1.0,1.0)[0]/quad(lambda x:f_cos(x)**2,-1.0,1.0)[0]
# error_2 = quad(lambda x:(c[0]*y_0(x)+c[1]*y_1(x)+c[2]*y_2(x)-f_cos(x))**2,-1.0,1.0)[0]/quad(lambda x:f_cos(x)**2,-1.0,1.0)[0]
# error_1 = quad(lambda x:(c[0]*y_0(x)+c[1]*y_1(x)-f_cos(x))**2,-1.0,1.0)[0]/quad(lambda x:f_cos(x)**2,-1.0,1.0)[0]
# error_0 = quad(lambda x:(c[0]*y_0(x)+0*x-f_cos(x))**2,-1.0,1.0)[0]/quad(lambda x:f_cos(x)**2,-1.0,1.0)[0]
# print(error_0)
#
# print(error_1)
# print(error_2)
# print(error_3)
# print(error_4)
# print(error_5)
#
#
# plt.semilogy([0, 1, 2,3,4,5],[error_0,error_1,error_2,error_3,error_4,error_5])
# plt.xlabel("k")
# plt.ylabel("log relative  error")
# plt.show()
#
#
#

