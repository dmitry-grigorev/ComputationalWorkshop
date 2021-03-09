from math import cos, pi
from mft1_2 import get_values, get_ERs, get_Ns, get_Pks


def equidistant_interpolation_Newton(f, pivot, a, b, degree, h):
    values = get_values(f, a, b, h)

    ERs = get_ERs(values, degree)

    N = get_Ns(a, b, degree, h, pivot)

    P, ER0 = get_Pks(a, b, degree, h, pivot, values, N, ERs)
    return P[len(P) - 1]


def non_equidistant_interpolation_Newton(f, nodes, pivot, degree):
    sorted_nodes = (sorted(nodes, key=lambda x: abs(x - pivot)))  # узлы в порядке близости к точке интерполяции
    DDs = get_DDs(f, sorted_nodes, degree)
    value = 0
    buf = 1
    for i in range(degree):
        value = value + DDs[i][0] * buf
        buf *= (pivot - sorted_nodes[i])
    return value


def get_DDs(f, nodes, degree):
    DDs = [[]] * degree
    DDs[0] = [f(nodes[k]) for k in range(degree)]
    for i in range(degree - 1):
        DDs[i + 1] = ([(DDs[i][j + 1] - DDs[i][j]) / (nodes[i + j + 1] - nodes[j]) for j in range(degree - 1 - i)])
    return DDs
