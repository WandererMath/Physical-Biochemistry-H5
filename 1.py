from math import comb
f_common=[0.9893,
    
    0.9963,
    0.9976,
    0.99988
    ]

f=[
    0.011,
    0.00365,
    0.00037,
    0.00015
]

n=[40,9,16,65]
atoms=["C13",  "N15", "O17", "H2"]

def p1(j):
    r=1
    r*=n[j]*f[j]*f_common[j]**(n[j]-1)
    for i, f0 in enumerate(f_common):
        if j!=i:
            r*=f0**n[i]
    return r

def p2(j, k):

    r=1
    if j==k:
        r*=comb(n[j], 2)*f[j]**2*f_common[j]**(n[j]-2)
        for i, f0 in enumerate(f_common):
            if j!=i:
                r*=f0**n[i]
        return r
    else:
        r*=n[j]*f[j]*n[k]*f[k]*f_common[j]**(n[j]-1)*f_common[k]**(n[k]-1)
        for i, f0 in enumerate(f_common):
            if j!=i and k!=i:
                r*=f0**n[i]
        return r

if __name__=="__main__":
    b=0
    for j, atom in enumerate(atoms):
        print(f"Fraction of +1 neutron form peptide (+1 neutron from {atom}) = {p1(j):.4f}")
        b+=p1(j)
    print(b)

    b=0
    tmp=[]
    for j, atom in enumerate(atoms):
        for k, atom2 in enumerate(atoms):
            if [j, k] not in tmp and [k, j] not in tmp:
                tmp.append([j, k])
                print(f"Fraction of +2 neutron form peptide (+1 neutron from {atom}, {atom2}) = {p2(j, k):.6f}")
                b+=p2(j, k)
    print(b)

