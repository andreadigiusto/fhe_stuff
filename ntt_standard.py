import numpy as np
from numpy.polynomial import Polynomial

#functions for ntt splits

def split6(coeff,root6,p):
    #root6 is a 6th root of 1, p is a prime modulus; this function computes the first split needed for our ntt decomposition;
    #uses the CRT split of x^n-x^(n/2)+1 = (x^(n/2) - root6) * (x^(n/2) - root6^5); observe that root6*root6^5=1 and root6+root6^5=1;

    coeff1, coeff2 = [], []
    x = len(coeff) // 2

    for i in range(x):
        j = i + x
        prod = root6 * coeff[j]
        coeff1.append((coeff[i] + prod) % p)
        coeff2.append((coeff[i] + coeff[j] - prod) % p)

    return coeff1, coeff2

def split2(coeff,root2,p):
    #root2 is a 2th root of r where coeff is the list of coefficients of a polynomial in Z_p[x]/x^k-r, p is a prime modulus;
    # this function computes one of the intermediate radix-2 splits needed for our ntt decomposition;
    #uses the CRT split of x^k = (x^(k/2) - root2) * (x^k - (-root2)); observe that root2*(-root2)=r;

    coeff1, coeff2 = [], []
    x = len(coeff) // 2

    for i in range(x):
        j = i + x
        prod = root2 * coeff[j]
        coeff1.append((coeff[i] + prod) % p)
        coeff2.append((coeff[i] - prod) % p)

    return coeff1, coeff2

"""
def invsplit6(coeff1,coeff2,root6,det):
    list1, list2 = [], []
    x = len(coeff1)

    for i in range(x):
        j = i + x
        prod = root6 * coeff1[i]
"""

#functions for polynomial multiplications in small degree rings

def multdeg3(a,b,r,p):
	# multiplication of two polynomials modulo X^3-r
	res=[]
	res.append((a[0]*b[0]+(a[1]*b[2]+a[2]*b[1])*r) %p )
	res.append((a[0]*b[1]+a[1]*b[0]+a[2]*b[2]*r) % p)
	res.append((a[0]*b[2]+a[2]*b[0]+a[1]*b[1]) % p)
	return res


#-------------------------------------------------------------------------------------------------------------------------------------------------

#tests

    #1. this code completely splits Z_13[x]/Phi_12(x) down to degree 0 factors using polynomial splits
q = 13
m = 12
root_tree = [[4],[2,6]]
Phi_12 = [1,0,-1,0,1]

def ntt12(pol):
    pol0, pol1 = split6(pol,root_tree[0][0],q)
    pol00, pol01 = split2(pol0,root_tree[1][0],q)
    pol10, pol11 = split2(pol1,root_tree[1][1],q)

    return pol00 + pol01 + pol10 + pol11

A = [1,2,3,4]
B = [5,6,7,8]

def ntt12mult(polA,polB):
    nttA = ntt12(polA)
    nttB = ntt12(polB)
    return np.multiply(nttA, nttB) % q

#computations in magma show that the product of the two polynomials has coefficients [3,3,4,8]

#print(ntt12mult(A,B) == ntt12([3,3,4,8]))


    #2. this code splits Z_37[x]/Phi_36(x) down to factors of degree 3 (modulo x^3-r) where multiplication is performed in a textbook way

#the ring cannot be split completely because there are no 9th roots of unity in Z_37*

q = 37
m = 36
root_tree = [[11],[14,8]]
rlist = [14,23,8,29]
Phi_36 = [1,0,0,0,0,0,-1,0,0,0,0,0,1]       #x^12-x^6+1

def ntt36inc(pol):
    pol0, pol1 = split6(pol,root_tree[0][0],q)
    pol00, pol01 = split2(pol0,root_tree[1][0],q)
    pol10, pol11 = split2(pol1,root_tree[1][1],q)

    return [pol00, pol01, pol10, pol11]

def ntt36incmult(polA,polB):
    nttA = ntt36inc(polA)
    nttB = ntt36inc(polB)

    for i in range(len(nttA)):
        nttA[i] = multdeg3(nttA[i], nttB[i], rlist[i],q)
    
    return nttA

A = [1,2,3,4,5,6,7,8,9,0,9,8]
B = [7,6,5,4,3,2,1,2,3,4,5,6]

#computations in magma show that the product of the two polynomials has coefficients [ 28, 10, 5, 16, 30, 5, 20, 29, 11, 9, 30, 14 ]
print(ntt36incmult(A,B) == ntt36inc([ 28, 10, 5, 16, 30, 5, 20, 29, 11, 9, 30, 14 ]))