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

if __name__=="__main__":
    b=0
    for j, atom in enumerate(atoms):
        print(f"Fraction of +1 neutron form peptide (+1 neutron from {atom}) = {p1(j):.4f}")
        b+=p1(j)
    print(b)
