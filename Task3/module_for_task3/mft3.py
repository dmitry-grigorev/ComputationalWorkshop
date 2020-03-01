from pandas import DataFrame


def get_RRs(values):
    x_s = []
    for node in values.keys():
        x_s += [node] * len(values[node])
    l = len(x_s)
    col1 = [values[x][0] for x in x_s]
    col2 = [values[x_s[i]][1] if x_s[i] == x_s[i + 1]
            else (values[x_s[i + 1]][0] - values[x_s[i]][0]) / (x_s[i + 1] - x_s[i])
            for i in range(l - 1)]

    col3 = [values[x_s[i]][2] / 2 if x_s[i] == x_s[i + 2]
            else (col2[(2 * i + 3) // 2] - col2[(2 * i + 1) // 2]) / (x_s[i + 2] - x_s[i])
            for i in range(l - 2)]
    RRs = [col1, col2, col3]
    for i in range(3, l):
        RRs.append([(RRs[i - 1][k + 1] - RRs[i - 1][k]) / (x_s[k + i] - x_s[k])
                    for k in range(l - i)])
    return RRs


def Hermite_interp(values):
    x_s = []
    for node in values.keys():
        x_s += [node] * len(values[node])
    l = len(x_s)
    RRs = get_RRs(values)
    column1 = [x_s[i // 2] if i % 2 == 0 else ' ' for i in range(2 * l - 1)]
    column2 = [values[x_s[i // 2]][0] if i % 2 == 0 else ' ' for i in range(2 * l - 1)]
    othercols = [column1, column2]
    for i in range(1, l):
        othercols.append([' '] * i + [RRs[i][k // 2] if k % 2 == 0
                                      else ' ' for k in range(2 * (l - i - 1) + 1)] + [' '] * i)

    table = DataFrame(othercols, columns=[' '] * (2 * l - 1),
                      index=['$x$', '$f(x)$'] + ['р.р. {0:} пор.'.format(i + 1)
                                                 for i in range(l - 1)]).T
    return table, [RRs[i][0] for i in range(l)]
