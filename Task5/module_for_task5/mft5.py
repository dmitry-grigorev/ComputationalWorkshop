from mft4 import non_equidistant_interpolation_Newton
from pandas import DataFrame


def derivate_experiment(derivative, f, h, pivot, nodes, df):
    values = []
    errors = []
    h_s = [h]
    values.append(derivative(f, h, pivot, nodes))
    errors.append(df(pivot) - values[0])
    h /= 2
    values.append(derivative(f, h, pivot, nodes))
    errors.append(df(pivot) - values[1])
    h_s.append(h)
    k = 1
    while errors[k] - errors[k - 1] < 0:
        h /= 2
        k += 1
        h_s.append(h)
        values.append(derivative(f, h, pivot, nodes))
        errors.append(df(pivot) - values[k])
    return DataFrame([h_s, values, errors], index=["$h$", "Численное значение производной", "Погрешность"])


def get_table(f, nodes, degree, h, df, d2f):
    values_in_nodes = [f(x) for x in nodes]
    real_val_der1 = [df(x) for x in nodes]
    real_val_der2 = [d2f(x) for x in nodes]
    num_val_der1_1 = [(values_in_nodes[i + 1] - values_in_nodes[i]) / h for i in range(degree)]
    num_val_der1_1 += [num_val_der1_1[degree - 1]]
    num_val_der1_2 = [(-3 * values_in_nodes[0] + 4 * values_in_nodes[1] - values_in_nodes[2]) / (2 * h)] + \
                     [(values_in_nodes[i + 1] - values_in_nodes[i - 1]) / (2 * h) for i in range(1, degree)] + \
                     [(3 * values_in_nodes[degree] - 4 * values_in_nodes[degree - 1] + values_in_nodes[degree - 2]) / (
                             2 * h)]
    errors1_1 = [real_val_der1[i] - num_val_der1_1[i] for i in range(degree + 1)]
    errors1_2 = [real_val_der1[i] - num_val_der1_2[i] for i in range(degree + 1)]
    num_val_der2 = [''] + [(values_in_nodes[i + 1] - 2 * values_in_nodes[i] + values_in_nodes[i - 1]) / (h ** 2) for i
                           in range(1, degree)] + ['']
    errors2 = [''] + [real_val_der2[i] - num_val_der2[i] for i in range(1, degree)] + ['']
    return DataFrame([nodes, values_in_nodes, real_val_der1, num_val_der1_1, errors1_1, num_val_der1_2, errors1_2,
                      real_val_der2, num_val_der2, errors2],
                     index=["$x$", "$f(x)$", "$f'(x)$", "$\widetilde {f'}(x), \ O(h)$",
                            "погр., $\ O(h)$", "$\widetilde{\widetilde {f'}}(x), \ O(h^2)$",
                            "погр., $ \ O(h^2)$", "$f''(x)$", "$\widetilde {f''}(x),\ O(h^2)$",
                            "погр., $ \ O(h^2)$"], columns=[''] * (degree + 1)).T
