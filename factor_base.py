#creates the factor base of N up to B

def eulersCriterion(n, p): #!!!!!!MAKE MORE EFFICIENT
    n = n % p
    for i in range(2, p, 1):
        if ((i * i) % p == n):
            return True
    return False

def factorBase(n, b):

    prime = [True for i in range(n+1)] 

    p = 2
    while(p * p <= n):
        if (prime[p] == True):
            
            for i in range(p * p, n + 1, p):
                prime[i] = False
        p += 1
    
    c = []
    for p in range(2, b):

        if prime[p] and eulersCriterion(n, p):
            c.append(p)
            
    #if 2 not in c: c.insert(0, 2)
    return c

print(factorBase(N, B))