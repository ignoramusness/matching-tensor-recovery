"""Independent, integer-only check of the four sharp rank certificates.

The Jacobian is built by evaluating F(W+e)-F(W), not by importing the
Jacobian constructor used to discover these certificates. Every coordinate
of F is affine in each scalar parameter separately, so these differences
are the exact partial derivatives. Bareiss elimination computes integer
determinants; no floating point or primality assumption is used.
"""
from itertools import combinations, product
import json
from pathlib import Path


def determinant(a):
    a=[list(row) for row in a]
    n=len(a)
    assert all(len(row)==n for row in a)
    if n==0:
        return 1
    previous, sign=1,1
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if a[i][k]),None)
        if pivot is None:
            return 0
        if pivot!=k:
            a[k],a[pivot]=a[pivot],a[k]
            sign=-sign
        p=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                num=a[i][j]*p-a[i][k]*a[k][j]
                assert num % previous == 0
                a[i][j]=num//previous
            a[i][k]=0
        previous=p
    return sign*a[-1][-1]


def tensor(w,words):
    return [w[0,1,a,b]*w[2,3,c,d]
            +w[0,2,a,c]*w[1,3,b,d]
            +w[0,3,a,d]*w[1,2,b,c]
            for a,b,c,d in words]


def check(record):
    dims=record['dimensions']
    columns=[(i,j,a,b) for i,j in combinations(range(4),2)
             for a,b in product(range(dims[i]),range(dims[j]))]
    words=list(product(*(range(d) for d in dims)))
    w=dict.fromkeys(columns,0)
    for i,j,a,b,value in record['nonzero_entries']:
        w[i,j,a,b]=value
    f=tensor(w,words)
    derivatives=[]
    for key in columns:
        shifted=dict(w)
        shifted[key]+=1
        derivatives.append([x-y for x,y in zip(tensor(shifted,words),f)])
    matrix=list(map(list,zip(*derivatives)))
    rows,cols=record['minor_rows'],record['minor_columns']
    assert len(rows)==len(cols)==record['rank']
    det=determinant([[matrix[i][j] for j in cols] for i in rows])
    assert det != 0
    assert det % record['prime']==record['minor_determinant_mod_prime']

    kernel=[]
    for r in range(3):
        lam=[int(i==r)-int(i==3) for i in range(4)]
        kernel.append([(lam[i]+lam[j])*w[key]
                       for key in columns for i,j in [key[:2]]])
    for root in range(4):
        if dims[root]!=2:
            continue
        others=[i for i in range(4) if i!=root]
        a,b,c=others
        def incident(vertex,label,coord):
            return (w[root,vertex,label,coord] if root<vertex
                    else w[vertex,root,coord,label])
        v=dict.fromkeys(columns,0)
        for i,j,s in ((b,c,1),(a,c,-1),(a,b,1)):
            for x,y in product(range(dims[i]),range(dims[j])):
                v[i,j,x,y]=s*(incident(i,0,x)*incident(j,1,y)
                              -incident(i,1,x)*incident(j,0,y))
        kernel.append([v[key] for key in columns])
    assert all(sum(x*y for x,y in zip(row,v))==0
               for row in matrix for v in kernel)
    kr,kc=record['kernel_minor_rows'],record['kernel_minor_columns']
    kd=determinant([[kernel[i][j] for j in kc] for i in kr])
    assert kd != 0
    assert kd % record['prime']==record['kernel_minor_determinant_mod_prime']
    return {'dimensions':dims,'jacobian_minor_size':len(rows),
            'jacobian_minor_integer_determinant':det,
            'kernel_minor_size':len(kr),'kernel_minor_integer_determinant':kd,
            'sorted_kernel_minor_integer_determinant':determinant(
                [[kernel[i][j] for j in sorted(kc)] for i in sorted(kr)]),
            'omitted_output_rows':sorted(set(range(len(words)))-set(rows)),
            'omitted_parameter_columns':sorted(set(range(len(columns)))-set(cols)),
            'sorted_row_minor_integer_determinant':determinant(
                [[matrix[i][j] for j in sorted(cols)] for i in sorted(rows)])}


if __name__=='__main__':
    folder=Path(__file__).resolve().parent
    records=json.loads((folder/'four_certificates.json').read_text(encoding='utf-8'))
    result=[check(r) for r in records]
    (folder/'four_exact_verification.json').write_text(
        json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
