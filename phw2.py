import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.linalg as la
from copy import copy
import matplotlib.cm as cm
import matplotlib

from scipy.linalg import lu as lu1
from mpl_toolkits.mplot3d import Axes3D



def e_generator(i,j,size):
    length = (size-1)**2
    e = np.zeros(length)
    e[(i-1)*(size-1)+j-1] = 1
    return e

def A_generator(N):
    Q = np.zeros([(N-1)**2, (N-1)**2])

    for i in range((N-1)**2):
        Q[i][i] = 4.0
        if i % (N-1) != 0:
            Q[i][i-1] = -1.0;
        if i % (N-1) != N-2:
            Q[i][i+1] = -1.0;
        if i // (N-1) != 0.0:
            if(i == 1):
                print(i / (N-1))
            Q[i][(i-(N-1))] = -1.0;
        if i // (N-1) < N-2:
            Q[i][(i+(N-1))] = -1.0;


    return Q
    '''notice that Q is the matrix. Q graph is for plt.
    Black: value 4
    Gray: value -1
    White: value 0
    '''

def lu(A):
    '''actually this is p^-1'''
    n = A.shape[0]
    U = A.copy()
    L = np.eye(n, dtype=np.double)
    P = np.eye(n, dtype=np.double)
    #index of largest abs value

    for i in range(n-1):
        #find the largest pivot
        max = abs(U[i,i])
        index = i
        for k in range(i,n):
            if abs(U[k,i])>max:
                max = abs(U[k,i])
                index = k

        #permutation
        if index != i:
            #swap in U, notice that swap in python is weird
            U[[index,i],i:n] = U[[i,index],i:n]
            P[[index,i]]=P[[i,index]]
            L[[index,i],:i] = L[[i,index],:i]

        for j in range(i+1,n):
            L[j,i] = U[j,i] / U[i,i]
            U[j,i:] -= L[j,i]*U[i,i:]
            U[j,i] = 0
    return P,L,U

def plot_0_white(A):
    plt.imshow(A,cmap = 'PiYG',norm = matplotlib.colors.Normalize(vmin = -1,vmax = 1),interpolation = 'none')
    plt.show()

def find_nonzero_entries(A):
    n = A.shape[0]
    count = 0
    for i in range(n):
        for(j) in range(n):
            if A[i][j] != 0:
                count = count + 1
    return count

def find_lower_bandwidth(A):
    n = A.shape[0]
    lower_bandwidth = 0
    flag = -1
    for i in range(1,n):
        for j in range(0, n-i):
            if A[i+j][j] != 0:
                flag = 1
                break;
        if flag == 1:
            lower_bandwidth = lower_bandwidth + 1
            flag = -1
        else:
            return lower_bandwidth;
    return lower_bandwidth+1;
'''L is matrix, b is array, return value x is also array'''
def fsolve(L,b):
    n = L.shape[0]
    b[0] = b[0]/L[0,0]
    for i in range(1,n):
        for j in range(0, i):
            b[i] = b[i] - L[i,j]*b[j]

        b[i] = b[i]/L[i,i]
    return b

def bsolve(U,b):
    n = U.shape[0]
    b[n-1] = b[n-1]/U[n-1,n-1]
    for i in range(n-2,-1,-1):
        for j in range(n-1,i,-1):
            b[i] = b[i] - U[i,j]*b[j]
        b[i] = b[i]/U[i,i]
    return b

'''sol on question 1'''


Q = A_generator(10)

plt.imshow(Q,cmap = 'PiYG',norm = matplotlib.colors.Normalize(vmin = -1,vmax = 1),interpolation = 'none')
plt.show()

'''test on question 2, lu(A) is my func'''
'''
A = np.array([[2, 5, 8, 7], [5, 2, 2, 8], [7, 5, 6, 6], [5, 4, 4, 8]])
L,U,P = lu(A)
print(P)
print(L)
print(U)
'''

'''question3 uncomment this to see all pics'''
#
Q_10 = A_generator(10)
P1,L1,U1 = lu(Q_10)
num1 = find_nonzero_entries(L1)
lb1 = find_lower_bandwidth(L1)
print(num1)
print(lb1)

Q_20 = A_generator(20)
P2,L2,U2 = lu(Q_20)
num2 = find_nonzero_entries(L2)
lb2 = find_lower_bandwidth(L2)
print(num2)
print(lb2)

Q_30 = A_generator(30)
P3,L3,U3= lu(Q_30)
num3 = find_nonzero_entries(L3)
lb3 = find_lower_bandwidth(L3)


Q_40 = A_generator(40)
P4,L4,U4= lu(Q_40)
num4 = find_nonzero_entries(L4)
lb4= find_lower_bandwidth(L4)
print(num4)
print(lb4)

Q_50 = A_generator(50)
P5,L5,U5 = lu(Q_50)
num5 = find_nonzero_entries(L5)
lb5 = find_lower_bandwidth(L5)
print(num5)
print(lb5)
x = np.array([num1,num2,num3,num4,num5])
y = np.array([10,20,30,40,50])
plt.title("nonzeros over N")
plt.scatter(y, x, color="black")
plt.show()
x = np.array(([lb1,lb2,lb3,lb4,lb5]))
plt.title("bandwidth over N")
plt.scatter(y, x, color="black")
plt.show()
'''sanity test for q3, IGNORE FOLLOWING THINGS  '''
# print(find_nonzero_entries(A))
# print(find_lower_bandwidth(A))
'''question4'''
'''sanity check for bsolve, fsolve, egenerator, please ignore'''
i = 2
j = 2
size = 4
matrix_e = e_generator(i,j,size)
print(matrix_e)
a = np.array([[1, 3, 5, 9], [0, 1, 0, 11], [0, 0, 2, 9], [0, 0, 0, 4]])
b = np.array([4.0, 2.0, 4.0, 2.0])
sol = bsolve(a,b)
print(sol)
print(a.dot(sol))
size = 20
Q_50 = A_generator(size)
e = e_generator(5,5,size)
P,L,U  = lu(Q_50)
e = e.dot(P)
ux = fsolve(L,e)
x = bsolve(U,ux)
grid = np.zeros([(size-1),(size-1)])
ind = 0
for i in range(0,(size-1)):
    for j in range(0,(size-1)):
        grid[i,j] = x[ind]
        ind += 1
nx, ny = size-1, size-1
x = range(nx)
y = range(ny)

hf = plt.figure()
ha = hf.add_subplot(111, projection='3d')
X, Y = np.meshgrid(x, y)  # `plot_surface` expects `x` and `y` data to be 2D
ha.plot_surface(X, Y, grid)
plt.show()
