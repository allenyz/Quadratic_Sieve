import numpy as np
from math import exp, sqrt, log
import time

from factor_base import *
from tonelli_shanks import *
from helpers import *

def smoothness_bound(n):
    B = 1.5*exp(sqrt(log(n) * log(log(n))) * 0.5)
    return int(B)

def find_smooth(n, b, fb, i):
    root = math.ceil(n**.5)

    a = []
    smoothx = [] #smooth number x's
    i_num = 0

    while(len(smoothx) < len(fb)):
        #print(i_num)
    
        #making relationships (x+k)^2 mod n
        a = interval(n, i_num*i + i, i_num*i)
        #print(a)
        
        for p in range(len(fb)): #prime
            for j in tonelli(n, fb[p]): #starting sequence number
                for k in range((j - i_num*i) % fb[p], len(a), fb[p]):
                    while a[k] % fb[p] == 0: #account for prime powers
                        a[k] //= fb[p]
        
    
        for o in range(len(a)):
            if a[o] == 1: smoothx.append(i_num*i + o)

        i_num+=1

    #print(smoothx)
    
    #for i in smoothx:
    #    print((root + i)**2 - n)
    #print(len(smoothx))
    #quit()

    
    return smoothx

def create_matrix(n, fb, smooth_x):

    parity_matrix = []
    
    for s in smooth_x:
        smooth_number = x_to_number(n, s)
        parity_row = []
        for p in fb:
            parity_row.append(0)
            while (smooth_number % p == 0):
                parity_row[-1]^=1 #flips parity
                smooth_number //= p
        parity_matrix.append(parity_row)

    return parity_matrix


def assert_larger_interval(fb, smooth):
    if (len(smooth) < len(fb)):
        print("Not Enough Smooth Numbers - Increase Interval Size")
        quit()



'''A Fast Algorithm for
        Gaussian Elimination over GF (2) and its Implementation on the
        GAPP.' Journal of Parallel and Distributed Computing 13.1
        (1991): 118-122.'''

def gaussian_elim(mat_in):

    mat = np.array(mat_in)
    pRows = np.zeros(mat.shape[0])

    for i in range(mat.shape[1]):
        
        if (sum(mat[:, i]) != 0):
            
            p = list(mat[:, i]).index(1)
            pRows[p] = 1
            
            for j in range(mat.shape[1]):
                
                if (j == i):
                    continue
                    
                if (mat[p][j] == 1):
                    mat[:, j] = (mat[:, j] + mat[:, i]) % 2

    return mat, pRows #as np.arrays

#t1 = [[1, 1, 0, 0], [1, 1, 0, 1], [0, 1, 1, 1], [0, 0, 1, 0], [0, 0, 0, 1]]
#test = np.array(t1)

def get_relationship_rows(mat, pRows):
    potentialRelationships = []
    rRows = []
    for i in range(len(pRows)):
        if(pRows[i] == 0):
            if (sum(mat[i]) != 0): #then potential relationship row
                potentialRelationships.append(mat[i])
                rRows.append(i)
    print(f"Got {len(potentialRelationships)} potential relations")
    return potentialRelationships, rRows


def get_dependent_rows(potentialRelationships, mat, index, pRows, rRows):
    potentialRelation = potentialRelationships[index]
    dependentRows = []
    dependentRows.append(rRows[index])
    for i in range(len(potentialRelation)):
        if (potentialRelation[i] == 1):
            for j in range(len(mat)):
                if (mat[j][i] == 1 and pRows[j] == 1):
                    dependentRows.append(j)
    return set(dependentRows)


def find_factors(n, fb, dependents, interval, smooth_x):
    root = math.ceil(n**.5)

    right_hand = 1
    left_hand = 1

    for i in dependents:
        left_hand *= ((root + smooth_x[i])**2-n)
        right_hand = (right_hand * (root + smooth_x[i]))
    
    left_hand = isqrt(left_hand)

    factor1 = gcd((left_hand - right_hand) % n, n)
    factor2 = n / factor1
    
    return factor1, factor2


def quadratic_sieve(n):
    B = smoothness_bound(n)
    #print(B)
    #B = 8000
    print(f"Looking for Factor Base of Size: {B}")
    I = 1000000
    iterator_num = 0

    fb = factorBase(n, B)
    print(f"Factor Base of Size {len(fb)}")

    sieve_interval = interval(n, I) #pre-smoothed number sequence
    #print(sieve_interval)

    smooth = find_smooth(n, B, fb, I) #indices of the smooth numbers from number sequence
    print(f"{len(smooth)} Smooth Relations")
    #print(len(fb))
    #print(len(smooth))
    #assert_larger_interval(fb, smooth)

   
    parity_matrix = create_matrix(n, fb, smooth) #large matrices with rows (smooth numbers) and columns (factor base)

    reduced_matrix, pRows = gaussian_elim(parity_matrix) #parity matrix after Gaussian elim, marked pivot rows

    potentialRelations, rRows = get_relationship_rows(reduced_matrix, pRows) #matrix of nonzero, nonpivot rows

    for i in range(len(potentialRelations)):
        dependentRows = get_dependent_rows(potentialRelations, reduced_matrix, i, pRows, rRows) #array of indices of dependent rows in reduced matrix
        factor1, factor2 = find_factors(n, fb, dependentRows, sieve_interval, smooth)
        if (factor1 != 1 and factor2 != 1):
            return factor1, factor2
    
    print("FACTORS NOT FOUND")
    return 1, n


if __name__ == "__main__":
    '''
    tests = [16921456439215439701, 
             46839566299936919234246726809,
             6172835808641975203638304919691358469663,
             3744843080529615909019181510330554205500926021947]
    '''
    tests = [6172835808641975203638304919691358469663]
    
    for test in tests:
        start_time = time.time()
        print(quadratic_sieve(test))
        print("Time: %s seconds" % (time.time() - start_time))