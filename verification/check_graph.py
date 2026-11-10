"""Independent graph-support audit, rebuilt from the definition of the map.

No manuscript derivative formula or archived source code is imported.
The prime-modular checks supplement the analytic audit, not replace it.
"""
from collections import defaultdict
from itertools import combinations, product
from pathlib import Path
import json
import random

PRIME = 1000003

def matchings(vertices):
    if not vertices:
        yield ()
        return
    i = vertices[0]
    for k, j in enumerate(vertices[1:], 1):
        for tail in matchings(vertices[1:k] + vertices[k+1:]):
            yield ((i,j),) + tail

def rank_sparse(columns):
    pivots = {}
    for column in columns:
        col = {key: value % PRIME for key, value in column.items() if value % PRIME}
        while col:
            p = min(col)
            if p not in pivots:
                inverse = pow(col[p], -1, PRIME)
                pivots[p] = {key: value * inverse % PRIME for key, value in col.items()}
                break
            factor = col[p]
            for key, value in pivots[p].items():
                new = (col.get(key, 0) - factor * value) % PRIME
                if new:
                    col[key] = new
                else:
                    col.pop(key, None)
    return len(pivots)

def matrix_columns(dims, pm, variable_edges):
    """Expand each matching product with one basis-coordinate variation.

    Each undifferentiated edge block is the three-dimensional identity.
    The selected coefficient is exactly the coefficient of t in the
    matching-product expression after replacing the selected block by I+tE.
    """
    n = len(dims)
    columns = []
    for edge in variable_edges:
        i, j = edge
        containing = [m for m in pm if edge in m]
        for a in range(dims[i]):
            for b in range(dims[j]):
                col = defaultdict(int)
                for matching in containing:
                    rest = [e for e in matching if e != edge]
                    for labels in product(range(3), repeat=len(rest)):
                        word = [None] * n
                        word[i], word[j] = a, b
                        for (u,v), label in zip(rest, labels):
                            word[u] = word[v] = label
                        col[tuple(word)] += 1
                columns.append(col)
    return columns

def check(dims, edges, variables=None, label=''):
    edges = set(edges)
    pm = [m for m in matchings(tuple(range(len(dims)))) if set(m) <= edges]
    active = set().union(*(set(m) for m in pm)) if pm else set()
    if variables is None:
        variables = edges
    variables = sorted(set(variables))
    effective = sorted(set(variables) & active)
    incidence_columns = [{r:1 for r,m in enumerate(pm) if e in m} for e in effective]
    incidence_rank = rank_sparse(incidence_columns)
    expected = sum(dims[i]*dims[j] for i,j in effective)-len(effective)+incidence_rank
    obtained = rank_sparse(matrix_columns(dims, pm, variables))
    assert obtained == expected, (label, dims, edges, variables, obtained, expected)
    return dict(label=label, dimensions=dims, edges=sorted(edges),
                variable_edges=variables, active_edges=sorted(active),
                perfect_matchings=len(pm), incidence_rank=incidence_rank,
                rank_mod=obtained, expected=expected)

def main():
    checks = []
    edges4 = list(combinations(range(4),2))
    for mask in range(1 << len(edges4)):
        edges = [e for k,e in enumerate(edges4) if mask >> k & 1]
        checks.append(check((3,3,3,3), edges, label=f'all-four-supports-{mask}'))
    edges6 = list(combinations(range(6),2))
    pm6 = list(matchings(tuple(range(6))))
    rng = random.Random(19071)
    supports = [set(edges6), {(i,(i+1)%6) if i<(i+1)%6 else ((i+1)%6,i) for i in range(6)},
                {(i,j) for i in range(3) for j in range(3,6)}, set(pm6[0]),
                {(0,1),(0,2),(1,2),(3,4),(3,5),(4,5)},
                set(pm6[0]) | {(0,2),(2,4)}]
    for count in [2,3,4,6,9,12]:
        supports.append(set().union(*(set(m) for m in rng.sample(pm6,count))))
    for k, edges in enumerate(supports):
        checks.append(check((3,3,3,3,3,3),edges,label=f'six-support-{k}'))
        variables = [e for e in sorted(edges) if rng.randrange(2)]
        checks.append(check((3,3,3,3,3,3),edges,variables,label=f'six-selected-{k}'))
    checks.append(check((3,4,3,5),edges4,label='unequal-four-complete'))
    checks.append(check((3,4,3,4,3,3),supports[1],label='unequal-six-cycle'))
    result = dict(prime=PRIME,definition='coefficient of t in matching product at padded I3',
                  count=len(checks),checks=checks)
    output = Path(__file__).with_name('graph-check-results.json')
    output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=len(checks),all_passed=True,output=output.name)))

if __name__ == '__main__':
    main()
