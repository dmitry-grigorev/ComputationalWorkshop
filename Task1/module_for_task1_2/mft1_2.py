from mft1 import find_max_abs
import pandas as pd


def get_values(f, a, b, h):  # Вычисляет значения функции в узлах интерполирования
    n = int((b - a) / h) + 1
    return [f(a + k * h) for k in range(n)]


def buildTable1(a, b, h, values):  # Строит таблицу точка-значение
    n = int((b - a) / h) + 1
    nodes = [(a + k * h) for k in range(n)]
    return pd.DataFrame([nodes, values], index=["x", "f(x)"], columns=[" "] * n).T


def get_ERs(values, degree):  # Вычисляет конечные разности
    le = len(values)
    ERs = [0]
    ERs[0] = [values[k + 1] - values[k] for k in range(le - 1)]
    for i in range(1, degree):
        ERs.append([ERs[i - 1][k + 1] - ERs[i - 1][k] for k in range(le - 1 - i)])
    return ERs


def buildTable2(a, b, degree, h, values, ERs):  # Строит таблицу конечных разностей
    n = int((b - a) / h) + 1
    nodes = [a + (k // 2) * h if k % 2 == 0 else ' ' for k in range(2 * n - 1)]
    values_table = [values[k // 2] if k % 2 == 0 else ' ' for k in range(2 * n - 1)]
    ER = []
    for i in range(degree):
        ER.append([' '] * i + [round(ERs[i][k // 2], 5)
                               if k % 2 == 1 else ' ' for k in range(2 * (n - i) - 1)] + [' '] * i)
    return pd.DataFrame([nodes, values_table] + ER,
                        index=["x", "y", '$ \Delta y $'] + ['$ \Delta^{0:} y $'.format(i+1) for i in range(1, degree)],
                        columns=[" "] * (2 * n - 1)).T


def get_Ns(a, b, degree, h, pivot):  # Вычисляет значения Nk(t)
    N = [1] + [0] * (degree + 1)
    if pivot <= a + h / 2:
        t = (pivot - a) / h
        for k in range(1, degree + 2):
            N[k] = N[k - 1] * (t - k + 1) / k
    elif pivot >= b - h / 2:
        t = (pivot - b) / h
        for k in range(1, degree + 2):
            N[k] = N[k - 1] * (t + k - 1) / k
    else:
        mid = a + h * int((b - a) / h) / 2
        if mid < pivot <= mid + h / 2:
            t = (pivot - mid) / h
            for k in range(1, degree + 2):
                N[k] = N[k - 1] * (t - (-1) ** k * k // 2) / k
    return N


def get_Pks(a, b, degree, h, pivot, values, N, ERs):  # Вычисляет значения Pk(pivot)
    Pks = [0] * (degree + 1)
    ER0 = [0] * (degree + 1)
    if pivot <= a + h / 2:
        ER0 = [ERs[k][0] for k in range(degree)]
        Pks[0] = values[0]
    else:
        n = len(values)
        if pivot >= b - h / 2:
            ER0 = [ERs[k][n - 2 - k] for k in range(degree)]
            Pks[0] = values[n - 1]
        else:
            mid = a + h * int((b - a) / h) / 2
            if mid < pivot <= mid + h / 2:
                ER0 = [ERs[k][(n - 1 - k) // 2] for k in range(degree)]
                Pks[0] = values[n // 2]
    for k in range(0, degree):
        Pks[k + 1] = Pks[k] + ER0[k] * N[k + 1]
    return Pks, [values[0]] + ER0


def buildTable3(f, derivatives, pivot, a, b, degree, h, ER0, N, P):  # Строит таблицу 1 по образцу
    string1 = ER0

    string2 = N[:-1]

    string3 = [ER0[k] * N[k] for k in range(degree + 1)]

    string4 = P

    string5 = [f(pivot) - P[k] for k in range(degree + 1)]

    maximums = [0] * (degree + 1)
    if pivot <= a + h / 2:
        maximums[0] = abs(derivatives[0](a))
        for k in range(1, min(degree + 1, len(derivatives))):
            maximums[k] = find_max_abs(derivatives[k], [a, a + k * h])

    else:
        if pivot >= b - h / 2:
            maximums[0] = abs(derivatives[0](b))
            for k in range(1, min(degree + 1, len(derivatives))):
                maximums[k] = find_max_abs(derivatives[k], [b - k * h, b])

        else:
            mid = a + h * int((b - a) / h) / 2
            if mid < pivot <= mid + h / 2:
                maximums[0] = abs(derivatives[0](mid))
                for k in range(1, min(degree + 1, len(derivatives))):
                    maximums[k] = find_max_abs(derivatives[k],
                                               [mid - (k // 2) * h, mid + ((k + 1) // 2) * h])
    string6 = [0] * (degree + 1)
    buf = h
    for k in range(0, degree + 1):
        string6[k] = maximums[k] * abs(N[k + 1]) * buf
        buf *= h

    dataframe = pd.DataFrame([string1, string2, string3, string4, string5, string6],
                             index=['$ \Delta^{k} y_{0} $', "$N_{k}(t)$",
                                    "$ N_{k}(t) \cdot \Delta^{k} y_{0} $",
                                    "$P_{0:}({1:})$".format("k", pivot),
                                    "$f({0:})-P_{1:}({2:})$".format(pivot, "k", pivot),
                                    "$|R_{0:}({1:})|\le$".format("k", pivot)],
                             columns=[str(k) for k in range(degree + 1)])
    dataframe.columns.name = "k"
    return dataframe


def equidistant_interpolation_tables(f, pivot, a, b, degree, h, derivatives):
    values = get_values(f, a, b, h)

    dataframe1 = buildTable1(a, b, h, values)

    ERs = get_ERs(values, degree)

    dataframe2 = buildTable2(a, b, degree, h, values, ERs)

    N = get_Ns(a, b, degree, h, pivot)

    P, ER0 = get_Pks(a, b, degree, h, pivot, values, N, ERs)

    dataframe3 = buildTable3(f, derivatives, pivot, a, b, degree, h, ER0, N, P)
    return dataframe1, dataframe2, dataframe3
