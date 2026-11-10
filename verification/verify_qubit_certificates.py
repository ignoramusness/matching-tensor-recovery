"""Exact integer certificates for the all-order qubit differential proof.

No third-party packages, random choices, data files, or modular assumptions.
Vertices in this program are numbered 1 through n.
For diagonal blocks diag(a,b), a labels 0 and b labels 1.
The directed variable y_ij has label 1 at i and label 0 at j.
"""
from functools import lru_cache
from itertools import combinations


def determinant(matrix):
    """Fraction-free Bareiss elimination, with checked exact divisions."""
    a = [list(row) for row in matrix]
    n = len(a)
    assert all(len(row) == n for row in a)
    sign, previous = 1, 1
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        value = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = value * a[i][j] - a[i][k] * a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
            a[i][k] = 0
        previous = value
    return sign * a[-1][-1]


def hafnian_function(edges):
    """Scalar hafnian for a binary, zero-diagonal symmetric graph matrix."""
    @lru_cache(None)
    def h(vertices):
        if not vertices:
            return 1
        assert len(vertices) % 2 == 0
        i = vertices[0]
        return sum(
            int(tuple(sorted((i, j))) in edges)
            * h(tuple(v for v in vertices if v not in (i, j)))
            for j in vertices[1:]
        )
    return h


def six_matrix():
    vertices = tuple(range(1, 7))
    a_edges = {(1,2),(1,4),(2,4),(2,6),(3,4),(3,5),(3,6),(4,6),(5,6)}
    b_edges = {(1,3),(1,5),(1,6),(2,4),(2,5),(2,6),(3,4),(4,5)}
    ha, hb = hafnian_function(a_edges), hafnian_function(b_edges)
    omitted = {(3,4,6),(3,5,6)}
    rows = [s for k in (1,3,5) for s in combinations(vertices, k)
            if s not in omitted]
    columns = [(i,j) for i in vertices for j in vertices if i != j]
    assert len(rows) == len(columns) == 30
    return [[
        hb(tuple(v for v in s if v != i))
        * ha(tuple(v for v in vertices if v not in s and v != j))
        if i in s and j not in s else 0
        for i,j in columns] for s in rows]


def eight_matrix():
    vertices = tuple(range(1, 9))
    b_edges = {(1,2),(1,4),(1,7),(1,8),(2,3),(2,6),(2,7),(2,8),
               (3,5),(3,6),(3,7),(3,8),(4,5),(4,7),(4,8),
               (5,6),(5,8),(6,8)}
    rows = list(combinations(vertices, 3))
    columns = [(i,j) for i in vertices for j in vertices if i != j]
    assert len(rows) == len(columns) == 56
    return [[int(tuple(v for v in s if v != i) in b_edges)
             if i in s and j in s else 0
             for i,j in columns] for s in rows]


if __name__ == '__main__':
    six = determinant(six_matrix())
    eight = determinant(eight_matrix())
    assert six == 128, six
    assert eight == -262144, eight
    print('Six-vertex 30 x 30 odd-block minor: determinant =', six)
    print('Eight-vertex 56 x 56 triangle matrix: determinant =', eight)
    print('Both exact nonzero integer certificates verified.')
