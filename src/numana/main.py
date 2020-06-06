from . import rootsolvers
import logging
import math

logging.basicConfig(level=logging.DEBUG)


def f(x):
    return math.atan(x)


def Df(x):
    return 1 / (1 + x * x)


def DDf(x):
    return (-2 * x) / (1 + x * x) ** 2



def main():
    a = rootsolvers.NR(f, Df, -1.392, 1E-5)
    print(f"f({a}) = {f(a)}")
    # error = 1E-5
    # start, stop = (-1,3)
    # i, i1 = (-1, 3)
    # a = rootsolvers.bisect(f, start, stop, error)
    # print(f"f({a}) = {f(a)}")
    # a = rootsolvers.RF(f, start, stop, error)
    # print(f"f({a}) = {f(a)}")
    # a = rootsolvers.fixed_point_iter(f, -1, error)
    # print(f"f({a}) = {f(a)}")
    # a = rootsolvers.secant(f, i, i1, error)
    # print(f"f({a}) = {f(a)}")
    # a = rootsolvers.NR(f, Df, i, error)
    # print(f"f({a}) = {f(a)}")

    # a = rootsolvers.mod_NR(f, Df, DDf, i, error)
    # print(f"f({a}) = {f(a)}")

