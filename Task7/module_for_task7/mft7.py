from math import sqrt
from pandas import DataFrame


def calculate_det(a, b, c, d):  # Вычисляет определитель матрицы 2x2
    return a * d - b * c


def get_dk(f11, f12, f21, f22, x, y):  # Вычисляет d^(k), dx^(k), dy^(k)
    return calculate_det(f11(x, y), f12(x, y), f21(x, y), f22(x, y))


def norm(x0, y0, x1, y1):  # Вычисляет норму вектора
    return sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2)


def dataframe_row(xk, yk, fk, gk, norm_val):  # строит строчку таблицы
    return DataFrame([[xk, yk, norm_val, fk, gk]], columns=["$x_k$", "$y_k$", "$||(x_k-x_{k-1}),(y_k-y_{k-1})||$",
                                                            "$f(x_k,y_k)$", "$g(x_k,y_k)$"])


def Newton_method(f, g, df, dg, x0, y0, imax, eps):  # итоговая функция
    dataframe = dataframe_row(x0, y0, f(x0, y0), g(x0, y0), "$-$")
    dataframe.columns.name = "$k$"

    x1 = x0
    y1 = y0
    x0 = x1 + 1
    y0 = y1 + 1  # Всё это сделано для того, чтобы войти в цикл и не сломать логику в нём
    i = 0
    norma = norm(x0, y0, x1, y1)
    while norma > eps and i < imax:
        x0 = x1
        y0 = y1
        dfk = get_dk(f, df['y'], g, dg['y'], x0, y0)
        dgk = get_dk(df['x'], f, dg['x'], g, x0, y0)
        dk = get_dk(df['x'], df['y'], dg['x'], dg['y'], x0, y0)
        x1 = x0 - dfk / dk
        y1 = y0 - dgk / dk
        norma = norm(x0, y0, x1, y1)
        i += 1
        dataframe = dataframe.append(dataframe_row(x1, y1, f(x1, y1), g(x1, y1), norma), ignore_index=True)
    return dataframe, [x1, y1]  # даже,когда максимальное число итераций наступит, вернётся xkmax,ykmax
