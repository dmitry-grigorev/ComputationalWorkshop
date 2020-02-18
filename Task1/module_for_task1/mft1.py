import pandas as pd
import scipy.optimize


def getRRs(nodes, f):  # Функция получения списков разделённых разностей
    RR1, RR2, RR3 = [], [], []
    l = len(nodes)
    for i in range(l - 1):  # Получаем набор раздел. разн. 1-ого пор.
        delta = nodes[i + 1] - nodes[i]
        RR1.append((f(nodes[i + 1]) - f(nodes[i])) / delta)
    for i in range(l - 2):  # Получаем набор раздел. разн. 2-ого пор.
        delta = nodes[i + 2] - nodes[i]
        RR2.append((RR1[index_RR1(i + 1)] - RR1[index_RR1(i)]) / delta)

    for i in range(l - 3):  # Получаем набор раздел. разн. 3-ого пор.
        delta = nodes[i + 3] - nodes[i]
        RR3.append((RR2[index_RR2(i + 1)] - RR2[index_RR2(i)]) / delta)
    return RR1, RR2, RR3


def index_RR1(i):  # Функции навигации по спискам р.р.
    return (2 * i + 1) // 2


def index_RR2(i):
    return (index_RR1(i) + index_RR1(i + 1)) // 2


def index_RR3(i):
    return (index_RR2(i) + index_RR2(i + 1)) // 2


def buildTable1(pivot, f, nodes, RR1, RR2, RR3):  # Строит таблицу 1 по образцу
    column1 = [nodes[i // 2] if i % 2 == 0 else ' ' for i in range(2 * len(nodes) - 1)] + [pivot]

    column2 = [round(f(nodes[i // 2]), 5) if i % 2 == 0 else ' ' for i in range(2 * len(nodes) - 1)] + ['?']

    column3 = [round(RR1[i // 2], 5) if i % 2 == 1 else ' '
               for i in range(2 * len(nodes) - 1)] + [' ']

    column4 = [' '] + [round(RR2[i // 2], 5) if i % 2 == 1 else ' '
                       for i in range(2 * len(nodes) - 3)] + [' '] * 2

    column5 = [' '] * 2 + [round(RR3[i // 2], 3) if i % 2 == 1 else ' '
                           for i in range(2 * len(nodes) - 5)] + [' '] * 3

    return (pd.DataFrame((column1, column2, column3, column4, column5),
                         index=['x', 'f(x)', 'р.р 1 п', 'р.р 2 п', 'р.р 3  п'],
                         columns=[' '] * (2 * len(nodes)))).T


def buildTable2(f, pivot, neighbors, polynom, maximums, errors):  # Строит таблицу 1 по образцу
    string1 = [i for i in range(4)]

    string2 = neighbors

    string3 = polynom

    string4 = [f(pivot) - polynom[i] for i in range(4)]

    string5 = maximums

    string6 = errors

    return pd.DataFrame([string1, string2, string3, string4, string5, string6],
                        index=["i", "Узлы интерполирования в порядке их использования",
                               "$P_{i}($" + str(pivot) + "$) -$ Значение многочлена в точке интерполирования",
                               "$f($" + str(pivot) + "$) {-} P_{i}($" + str(pivot) + "$) -$ Фактическая погрешность",
                               "$M_{i+1} -$ Оценка модуля произодных",
                               "$R_{i}($" + str(pivot) + "$) -$ Оценка погрешности"], columns=[" "] * 4)


def interpolating_polynom(f, pivot, neighbors):  # Вычисление многочлена
    RR1_4 = [0] * 3
    for i in range(3):
        RR1_4[i] = (f(neighbors[i + 1]) - f(neighbors[i])) / (neighbors[i + 1] - neighbors[i])

    RR2_4 = [0] * 2
    for i in range(2):
        RR2_4[i] = (RR1_4[i + 1] - RR1_4[i]) / (neighbors[i + 2] - neighbors[i])

    RR3_4 = [0] * 1
    for i in range(1):
        RR3_4[i] = (RR2_4[i + 1] - RR2_4[i]) / (neighbors[i + 3] - neighbors[i])

    values = [0] * 4
    values[0] = f(neighbors[0])

    buf = (pivot - neighbors[0])
    values[1] = values[0] + RR1_4[0] * buf

    buf *= (pivot - neighbors[1])
    values[2] = values[1] + RR2_4[0] * buf

    buf *= (pivot - neighbors[2])
    values[3] = values[2] + RR3_4[0] * buf
    return values


def discard4(p):  # Округляет число с плавающей точкой до 4-ёх знаков
    return int(p * (10 ** 4)) / (10 ** 4)


def find_max_abs(f, bounds):  # Ищет максимум модуля функции на отрезке с заданными границами
    a = bounds[0]
    b = bounds[1]
    res = (-1) * (scipy.optimize.minimize_scalar(lambda x: -abs(f(x)), method='bounded', bounds=(a, b))).fun
    return discard4(res)


def find_estimation_der(list_of_der, pivot, neighbors):  # Ищет оценку модуля производных функции
    points = [pivot]
    maximums = [0] * len(list_of_der)
    for i in range(len(list_of_der)):
        points.append(neighbors[i])
        points = sorted(points)
        a = points[0]
        b = points[len(points) - 1]
        maximums[i] = find_max_abs(list_of_der[i], [a, b])
    return maximums


def interpolating(f, nodes, pivot, derivatives):  # Итоговая функция, решающая задачу
    RR1, RR2, RR3 = getRRs(nodes, f)

    dataframe1 = buildTable1(pivot, f, nodes, RR1, RR2, RR3)

    neighbors = (sorted(nodes, key=lambda x: abs(x - pivot)))[0:4]  # четыре ближайшие точки

    maximums = find_estimation_der(derivatives, pivot, neighbors)

    buf = 1
    errors = [0] * 4
    for k in range(4):
        buf *= abs(pivot - neighbors[k]) / (k + 1)
        errors[k] = maximums[k] * buf

    polynom = interpolating_polynom(f, pivot, neighbors)

    dataframe2 = buildTable2(f, pivot, neighbors, polynom, maximums, errors)
    RRs = [RR1, RR2, RR3]
    return RRs, dataframe1, dataframe2
