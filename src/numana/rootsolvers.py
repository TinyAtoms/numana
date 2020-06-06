import math
import logging  # I'll use debug to provide the extra info for my assignments


# I'll probably need these somewhere along the line
# import numpy as np
# import sympy


# BRACKETED ROOTFINDERS


def bisect(function, start=-100, stop=98, epsilon=5E-5):
    '''
    using bisect. kinda like binary search but for functions and rootfinding
    Guaranteed to converge
    TODO; consider writing the function such that the debug info is only computed when debug is on
    '''
    initial_error = abs((stop-start)/(start + stop)) * 100
    iterations = math.ceil(math.log2(initial_error / epsilon))
    f_low = function(start)
    last = start  # for debug reasons
    for i in range(iterations):
        guess = 0.5 * (start + stop)
        f_guess = function(guess)
        error = abs(guess - last)  # debug reasons
        logging.debug(
            f"{i}:[{start}, {stop}] -> guess: {guess}, Ea={error}, epsilon_a={error/guess}")
        last = guess  # debug reasons
        condition = f_guess * f_low
        if condition > 0:
            start = guess
            f_low = function(guess)
        elif condition < 0:
            stop = guess
        else:
            logging.info(f"bisect took {i} iterations")
            return guess
    logging.info(f"bisect took {i} iterations")
    return guess


def RF(function, low=-100, high=101, epsilon=5E-5):
    '''
    This one uses a smarter heuristic to guess
    the method used is called 
    The Illinois algorithm, modified regula falsi, modified false position
    guaranteed to converge, and faster than bisect
    '''
    error = 100
    f_low = function(low)
    f_high = function(high)
    low_stuck = False
    high_stuck = False
    i = 0
    while error > epsilon:
        guess = (low * f_high - high * f_low) / (f_high - f_low)
        f_guess = function(guess)
        condition = f_guess * f_low
        error = abs((high-low)/(high + low)) * 100
        logging.debug(
            f"{i}: [{low}, {high}] -> guess: {guess}, Ea={error * guess}, epsilon_a={error}")
        if condition > 0:
            low = guess
            f_low = f_guess
            low_stuck = False
            if high_stuck:
                f_high /= 2
            high_stuck = not high_stuck  # if high, set to low etc
        elif condition < 0:
            high = guess
            f_high = f_guess
            high_stuck = False
            if low_stuck:
                f_low /= 2
            low_stuck = not low_stuck
        else:
            logging.info(f"modified RF took {i} iterations")
            return guess
        i += 1
    logging.info(f"modified RF took {i} iterations")
    return guess


# OPEN METHODS


def fixed_point_iter(function, last_guess=0, epsilon=5E-5):
    '''
    this has linear convergence. see book pg 144 for more info 
    TODO: handle cases where guess = 0 and you get divby0 error
    '''
    def mod_function(x):
        return function(x) + x

    guess = mod_function(last_guess)
    error = 100
    i = 1
    while error > epsilon:
        logging.debug(
            f"{i}: last_guess= {last_guess}, guess = {guess}, Ea = {error * guess}, epsilon={error}")
        error = abs((guess - last_guess) / guess) * 100
        last_guess = guess
        guess = mod_function(last_guess)
        i += 1
        if i > 200:
            break
    logging.info(f"fixed point iter took {i} iterations")
    return guess


def secant(function, guess, last_guess=0, epsilon=5E-5):
    '''
    Newton's method is not guaranteed to converge. 
    It also makes a number of assumptions, and can fail for a number of reasons.
    check https://en.wikipedia.org/wiki/Newton%27s_method#Failure_analysis
    # TODO: a more elegant way to handle div by 0
    '''
    f_last = function(last_guess)
    f_guess = function(guess)
    error = 100
    i = 0
    while error > epsilon:
        i += 1
        try:
            next_guess = guess - \
                (f_guess * (last_guess - guess))/(f_last - f_guess)
        except ZeroDivisionError:
            next_guess = guess - (f_guess * (last_guess - guess)) / 0.0000001
        try:
            error = abs((next_guess - guess) / next_guess) * 100
        except ZeroDivisionError:
            error = abs((next_guess - guess) / 0.0000001) * 100
        logging.debug(
            f"{i}: last={last_guess}, guess={guess}, next={next_guess}, Ea={error * guess},  epsilonA={error}")
        last_guess = guess
        guess = next_guess
        f_last = function(last_guess)
        f_guess = function(guess)

    logging.info(f"secant method took {i} iterations")
    return guess


def NR(function, dfunction, initial_guess=0, epsilon=5E-5):
    '''
    Newton's method is not guaranteed to converge. 
    It also makes a number of assumptions, and can fail for a number of reasons.
    check https://en.wikipedia.org/wiki/Newton%27s_method#Failure_analysis
    # TODO: a more elegant way to handle div by 0
    '''
    logging.debug(
        f"calling NR(f, df, initial={initial_guess}, error={epsilon})")
    # I could just use last_guess in the definition, but i prefer to see initial_guess when the IDE shows me arghelp
    last_guess = initial_guess
    error = 100
    i = 0
    while error > epsilon:  # I have the sinking suspicion that we want Et here not epsilon_a
        try:
            guess = last_guess - function(last_guess) / dfunction(last_guess)
        except ZeroDivisionError:
            guess = last_guess - function(last_guess) / epsilon

        try:
            error = abs((guess - last_guess) / guess * 100)
        except ZeroDivisionError:
            error = abs((guess - last_guess) / epsilon) * 100
        logging.debug(
            f"\t {i}: last = {last_guess}, guess={guess}, Ea={error * guess}, epsilonA={error}")
        last_guess = guess
        i += 1
        if i == 200:
            break
    logging.info(f"NR took {i} iterations")
    return guess


def mod_NR(function, dfunction, ddfunction,  initial_guess=0, epsilon=5E-5):
    '''
    Newton's method but modified
    '''
    last_guess = initial_guess
    f_guess = function(last_guess)
    df_guess = dfunction(last_guess)
    error = 1
    i = 0
    while error > epsilon:  # FIXME : really beginning to doubt the usage of ea> error
        guess = last_guess - \
            ((f_guess * df_guess) / ((df_guess*df_guess) -
                                     f_guess * ddfunction(last_guess)))
        error = abs((guess - last_guess) / guess * 100)
        logging.debug(
            f"{i}: last = {last_guess}, guess={guess}, Ea={error}, epsilonA={error * guess}")
        f_guess = function(last_guess)
        df_guess = dfunction(last_guess)
        last_guess = guess
        i += 1
        if i > 200:
            break

    logging.info(f"NR took {i} iterations")
    return guess
