#!/usr/bin/env python
import math

import numpy
import numpy as np
import pandas as pd
import sympy
import scipy.optimize as opt
from sympy import *
import matplotlib.pyplot as plt

#The first version cqn solve simple questions that satisfies(2/9)
# 1.no root without sign change like x^2 at x=0
# 2. no root exactly on the interval x=a or x = b

#version 1.1: edit a logic chunk to check whether there is root without sign change 2/10
#root_without_sign_change
# for x^2 in [-1,1], the sturm will find out that there is one
# There are tolerance problems and I noticed that there might be problem if 2 roots are close to each other!

#version 1.2:  This version focus on the lucky action that we fit the zero point exactly on a/b
#it helps to solve some extreme things

#version 1.3 I believe the codes will pass the script.

#version 1.4 Add hw2


# This is a  function to find root
# p: narray, coefficient from p0 to p_n x^n
# a,b  construct the interval
# return a list of roots
def findroots(p, a, b):
    p = np.array(p)

    global roots
    roots = []
    sturm_func_arr = sturm_func(p)
    find_one_sign_change(sturm_func_arr,a,b)
    root_set = set(roots)
    unique_roots = list(root_set)
    unique_roots.sort()
    return unique_roots

def findroots2(p, a, b):
    tol = 1e-13
    p = np.array(p)

    roots = np.roots(p[::-1]) # np.roots takes coefficients in reverse order
    roots = [z for z in roots if abs(np.imag(z)) < tol
             and a - tol < np.real(z) < b + tol]
    return np.array(sorted(roots))
def find_one_sign_change(sturm_func, a, b):
    #1.test how many sign change
    #2.if change = 1, send to brentq
    #3.if sign = 0, neglect
    #4. if sign > 1, bisection
    if a>b:
        return
    tol = 1e-13
    tol_change = 1e-6
    if abs(np.polyval(org_poly,a)) < tol:
        if len(roots) == 0 or abs(a-roots[-1]) > 10**(-6):
            roots.append(a)
        find_one_sign_change(sturm_func, a+tol_change,b)
        return
    if abs(np.polyval(org_poly,b)) < tol:
        if len(roots) == 0 or abs(b-roots[len(roots) - 1]) > 10**(-6):
            roots.append(b)
        find_one_sign_change(sturm_func, a,b-tol_change)
        return



    count_a = 0
    count_b = 0
    pre_a = np.polyval(sturm_func[0],a)
    pre_b = np.polyval(sturm_func[0],b)

    for poly in sturm_func:
        cur_a = np.polyval(poly,a)
        cur_b = np.polyval(poly,b)
        if cur_a * pre_a < 0:
            count_a = count_a + 1
        if cur_b * pre_b < 0:
            count_b = count_b + 1
        if cur_a != 0:
            pre_a = cur_a
        if cur_b != 0:
            pre_b = cur_b
    sign_changes = abs(count_b-count_a)
    if sign_changes == 0 :
        return;
    elif sign_changes == 1:
        return brentq_find_root(a,b)
    else:
        find_one_sign_change(sturm_func,a,(a+b)/2)
        find_one_sign_change(sturm_func,(a+b)/2,b)



# given n array, get the sturm func
def sturm_func(p):
    sChain  = []
    #from p_n to p0
    p = np.flip(p)
    #poly1 - f_n-2, poly2 - f_n-1
    poly1 = np.poly1d(p)
    global org_poly
    org_poly = poly1
    poly2 = np.polyder(poly1)
    cur = poly2
    sChain.append(poly1)
    sChain.append(poly2)
    while cur.order != 0:
        #sturm itself
        cur = poly2*((np.polydiv(poly1,poly2))[0]) - poly1
        poly1 = poly2
        poly2 = cur
        sChain.append(cur)
    return sChain

#solve question like y = x^2
def root_without_sign_change(a,b):

    tol = 1e-13 # a logical tolerance, if (y(x)-0 < tol, we consider it as a  root!)
    if abs(a-b) < tol:
        a =  0
    elif abs(np.polyval(org_poly,a)) <= tol:
        roots.append(a)
        return True
    elif abs(np.polyval(org_poly,b)) <= tol:
        roots.append(b)
        return True
    else:
        if not root_without_sign_change(a,(a+b)/2):
            root_without_sign_change((a+b)/2,b);

#given n array, and a interval with only 1 root, return root
def brentq_find_root(a,b):
    if np.polyval(org_poly,a)*np.polyval(org_poly,b) > 0:
        #use pure bisection here
        root = root_without_sign_change(a,b)
    else:
        root = opt.brentq(f, a, b)
        if len(roots) == 0 or abs(root-roots[-1]) > 10**(-6):
            roots.append(root)
    return root;

def f(x):
    return np.polyval(org_poly, x)



#the following part is question2
#let assume we

#r is array length 3, t is just t, d is array length3, abc are const
#x0->k, y0->q
def get_c_with_point(k, q, a, b, c):
    t = symbols('t')
    sth = expand(k**4+q**4+(4-t)**4+a*(k**2+q**2+(4-t)**2)**2+b*(k**2+q**2+(4-t)**2)+c)
    p_of_t = Poly(sth)
    co = p_of_t.all_coeffs()
    co.reverse()
    return co
    # cof.append(sth.coeff(1))
    # cof.append(sth.coeff(t))
    # cof.append(sth.coeff(t**2))
    # cof.append(sth.coeff(t**3))
    # cof.append(sth.coeff(t**1))
    # print(cof)
if __name__ == '__main__':
    #hi, use this part to look at the solution
    #first line is my method
    #second line is the np root, with 0.j terms should appear in my result
    #be careful, the sequence of 2 array should not be same. reserve it!(seriously
    #
    roots_true = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    tol = 1e-13
    p = np.array([0.0, 0.09375, -0.78125, 2.1875, -2.5, 1.0])

    print("homework1, question1")
    print("this is the solution for my method! ")
    a = -0
    b = 1
    root_sol = findroots(p,a,b)
    print(root_sol)
    assert (abs(root_sol - roots_true) < tol*roots_true + tol).all()
    print(org_poly)

    #settings for question2
    #the camera is 4*4
    #z = 4, x:[-4,4] y:[-4,4]
    #direction[x=0,y=0,z=-1]
    #dx = 0, dy = 0, dz = -1
    #number represents the pixels for row/column, 512 -> 512*512 picture
    number = 512
    a = 0.0
    b = 0.0
    c = -1.0
    #separate the interval by 512 sections
    diff = 8.0/number
    z_0 = 4.0
    dx = 0.0
    dy = 0.0
    dz = -1.0
    #pixels
    #rx = x_0, ry = y_0, rz = z_0 - t = 4 - t
    color = 1
    #grid - 2D, stores C
    grid = []
    #grid is 512*512
    for i in range(number):
        new = []
        for j in range(number):
            new.append(0)
        grid.append(new)
    #add r_x,r_y into array
    arr_x = []
    arr_y = []
    for i in range(number):
        x = i*diff-4.0
        arr_x.append(x)
    for i in range(number):
        y = i*diff-4.0
        arr_y.append(y)
    i = 0
    j = 0
    for i in range(number):
        x_0 = arr_x[i]
        for j in range(number):
            y_0 = arr_y[j]

            #co: arr with coefficent
            co = get_c_with_point(x_0,y_0,a,b,c)
            co = np.array(co)
            #Please test my findroots separately. I discuss with Professor and next line fail
            #to run using findroots function, while typing coefficient array separately will make
            #findroots works well.
            res = findroots2(co,-4,4)
            if len(res) == 0:
                z_val  = 0
            else:
                z_val = 4-res[0]
            #find gradient
            vector = []
            vector.append(4*x_0**3+4*x_0*a*(x_0**2+y_0**2+z_val**2)+2*b*x_0)
            vector.append(4*y_0**3+4*y_0*a*(x_0**2+y_0**2+z_val**2)+2*b*y_0)
            vector.append(4*z_val**3+4*z_val*a*(x_0**2+y_0**2+z_val**2)+2*b*z_val)
            v = np.array(vector)


            #unit vector
            unit_vector = v / (v**2).sum()**0.5

            if(len(res) == 0):
                grid[i][j] = 0
            elif(np.real(unit_vector[2]) <= 1  and np.real(unit_vector[2]) >= -1):
                #arccos(-n*d)  while dx = 0, dy = 0,dz = -1
                #arccos(-n*d) = arccos(n_z) = arccos(unit  vector[2])

                grid[i][j] = math.acos(-np.real(unit_vector[2]))
            print(grid[i][j])

    plt.imsave('test.jpg',grid)
