#creates the factor base of N up to B

def eulersCriterion(n, p):
    
    ls = pow(n, (p - 1) // 2, p) 
    return ls if ls == 1 else -1

def factorBase(n, b):

    prime = [True for i in range(b+1)] 

    p = 2
    while(p * p <= b):
        if (prime[p] == True):
            
            for i in range(p * p, b + 1, p):
                prime[i] = False
        p += 1
    
    c = []
    for p in range(2, b):

        if (prime[p] and eulersCriterion(n, p) == 1):
            c.append(p)
            
    if 2 not in c: c.insert(0, 2)
    return c
