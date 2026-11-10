"""Fresh check of the tangent matrix at rank-three identity blocks.

No archived code or stored matrix is imported.  Formula comes directly from
pairing the remaining equal colours after differentiating one matching edge.
Modular ranks are lower bounds in characteristic zero; the proof gives equality.
"""
from itertools import product, combinations
import json
from pathlib import Path

PRIME = 1000003

def odd_df(count):
    out = 1
    for x in range(count-1, 0, -2):
        out *= x
    return out

def rank_mod(rows):
    basis = {}
    for original in rows:
        row = [v % PRIME for v in original]
        while any(row):
            lead = next(i for i, v in enumerate(row) if v)
            if lead not in basis:
                inv = pow(row[lead], -1, PRIME)
                basis[lead] = [v*inv % PRIME for v in row]
                break
            scale = row[lead]
            pivot = basis[lead]
            row = [(v-scale*w) % PRIME for v,w in zip(row,pivot)]
    return len(basis)

def kernel_dimension(n, selected):
    omitted = set(combinations(range(n), 2))-set(selected)
    constraints = [[1]*n]
    for i,j in omitted:
        row = [0]*n
        row[i] = row[j] = 1
        constraints.append(row)
    return n-rank_mod(constraints)

def check(dims, selected=None, baseline_rank=3):
    n = len(dims)
    selected = list(combinations(range(n),2)) if selected is None else selected
    cols = [(i,j,a,b) for i,j in selected for a in range(dims[i]) for b in range(dims[j])]
    def rows():
        for word in product(*(range(d) for d in dims)):
            row=[]
            for i,j,a,b in cols:
                if word[i]!=a or word[j]!=b:
                    row.append(0); continue
                rem=[word[k] for k in range(n) if k not in (i,j)]
                counts=[rem.count(c) for c in range(baseline_rank)]
                if any(c>=baseline_rank for c in rem) or any(k%2 for k in counts):
                    row.append(0); continue
                value=1
                for k in counts: value *= odd_df(k)
                row.append(value)
            yield row
    actual=rank_mod(rows())
    expected=len(cols)-kernel_dimension(n,selected)
    assert actual==expected, (dims,selected,actual,expected)
    return dict(dims=dims,baseline_rank=baseline_rank,edges=selected,columns=len(cols),rank_mod=actual,expected=expected)

def main():
    out=[]
    for dims in ([3]*4,[3]*6,[4]*4,[3,3,3,4]):
        out.append(check(dims))
    n=6
    for r in range(1,n+1):
        edges=[(i,j) for i,j in combinations(range(n),2) if i<r or j<r]
        out.append(check([3]*n,edges))
    # Complement graphs: balanced bipartite, imbalanced bipartite, odd cycle,
    # and disconnected mixtures test the general edge-subset formula.
    all_edges=set(combinations(range(n),2))
    omitted_sets=[[(0,1),(2,3),(4,5)],[(0,1),(0,2),(0,3),(0,4),(0,5)],
                  [(0,1),(1,2),(0,2)],[(0,1),(1,2),(0,3),(1,3),(2,3)],[]]
    for omitted in omitted_sets:
        out.append(check([3]*n,sorted(all_edges-set(omitted))))
    record={'prime':PRIME,'purpose':'Independent modular lower bounds; analytic proof supplies universal upper bounds','checks':out}
    path=Path(__file__).with_name('independent-rank-results.json')
    path.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'checks':len(out),'all_passed':True,'ranks':[x['rank_mod'] for x in out]}))

if __name__=='__main__': main()
