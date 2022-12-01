import numpy as np

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
    #root6 is a 2th root of r where coeff is the list of coefficients of a polynomial in Z_p[x]/x^k-r, p is a prime modulus;
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
