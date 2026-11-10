"""Independent integer certificate audit, built directly from matchings.

No hafnian-cofactor formula is used to build the derivative matrices.
Python standard library only. All vertex labels in this file are zero-based.
"""
from itertools import combinations
from hashlib import sha256
import json


def matchings(vertices):
    if not vertices:
        yield ()
        return
    i = vertices[0]
    for pos in range(1, len(vertices)):
        j = vertices[pos]
        rest = vertices[1:pos] + vertices[pos + 1:]
        for suffix in matchings(rest):
            yield ((i, j),) + suffix


def derivative_row(n, ones, a, b):
    """Differentiate each matching monomial at the binary diagonal source."""
    cols = [(i, j) for i in range(n) for j in range(n) if i != j]
    col_index = {edge: k for k, edge in enumerate(cols)}
    row = [0] * len(cols)
    ones = set(ones)
    for matching in matchings(tuple(range(n))):
        for edge_pos, (u, v) in enumerate(matching):
            if (u in ones) == (v in ones):
                continue
            i, j = (u, v) if u in ones else (v, u)
            weight = 1
            for other_pos, (s, t) in enumerate(matching):
                if other_pos == edge_pos:
                    continue
                if (s in ones) != (t in ones):
                    weight = 0
                    break
                weight *= (b if s in ones else a).get((s, t), 0)
            row[col_index[i, j]] += weight
    return row


def determinant(matrix):
    """Integer Bareiss elimination, checking every purported exact division."""
    a = [row[:] for row in matrix]
    n = len(a)
    assert all(len(row) == n for row in a)
    previous_pivot = 1
    sign = 1
    for k in range(n - 1):
        pivot_row = next((r for r in range(k, n) if a[r][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                quotient, remainder = divmod(numerator, previous_pivot)
                assert remainder == 0
                a[i][j] = quotient
            a[i][k] = 0
        previous_pivot = pivot
    return sign * a[-1][-1]


def fingerprint(matrix):
    return sha256(json.dumps(matrix, separators=(',', ':')).encode()).hexdigest()


def check_eight():
    n = 8
    b_one_based = [(1,2),(1,4),(1,7),(1,8),(2,3),(2,6),(2,7),(2,8),
                   (3,5),(3,6),(3,7),(3,8),(4,5),(4,7),(4,8),(5,6),
                   (5,8),(6,8)]
    b = {(i-1,j-1): 1 for i,j in b_one_based}
    a = {edge: 1 for edge in combinations(range(n), 2)}
    singleton_rows = [derivative_row(n, (i,), a, b) for i in range(n)]
    assert all(x % 15 == 0 for row in singleton_rows for x in row)
    normalized = [[x//15 for x in row] for row in singleton_rows]
    triangle_matrix = []
    for triangle in combinations(range(n), 3):
        direct = derivative_row(n, triangle, a, b)
        assert all(x % 3 == 0 for x in direct)
        row = [-x//3 for x in direct]
        for vertex in triangle:
            opposite_edge = tuple(x for x in triangle if x != vertex)
            weight = b.get(opposite_edge, 0)
            row = [x + weight*y for x,y in zip(row, normalized[vertex])]
        # Independently check the expected sparse description and ordering.
        for col, (i,j) in enumerate(( (i,j) for i in range(n)
                                      for j in range(n) if i != j )):
            wanted = 0
            if i in triangle and j in triangle:
                wanted = b.get(tuple(x for x in triangle if x != i), 0)
            assert row[col] == wanted
        triangle_matrix.append(row)
    value = determinant(triangle_matrix)
    assert value == -262144, value
    return {'shape': [56,56], 'determinant': value,
            'matrix_sha256': fingerprint(triangle_matrix)}


def check_six():
    n = 6
    a = {e:1 for e in [(0,1),(0,3),(1,3),(1,5),(2,3),(2,4),(2,5),(3,5),(4,5)]}
    b = {e:1 for e in [(0,2),(0,4),(0,5),(1,3),(1,4),(1,5),(2,3),(3,4)]}
    omitted = {(2,3,5),(2,4,5)}
    rows = [S for size in (1,3,5) for S in combinations(range(n), size)
            if S not in omitted]
    matrix = [derivative_row(n, S, a, b) for S in rows]
    value = determinant(matrix)
    assert value == 128, value
    return {'shape': [30,30], 'determinant': value,
            'matrix_sha256': fingerprint(matrix), 'omitted_rows': sorted(omitted)}


if __name__ == '__main__':
    print(json.dumps({'eight_vertex_triangle': check_eight(),
                      'six_vertex_odd_minor': check_six()}, indent=2))
