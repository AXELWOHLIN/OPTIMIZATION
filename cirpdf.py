import numpy as np
from scipy.stats import ncx2


def cirpdf(y, t, a, b, sigma):
    """
    Generate a vectorised CIR distribuion from noncentral chi-square distribution.

        dX = a(b - X)dt + sigma*sqrt(X) dW

    The CIR pdf is related to the Chi distribution with 2q+2 degrees of freedom and non-centrality parameter

    Args:
        y (numppy.array): a 1-d vector of evaluation points: [X_0, X_1, ..., X_n ]
        t (numpy.array): a 1-d vector of times:              [t_0, t_1, ..., t_n]
        a (float): parameter
        b (float): parameter
        sigma (float): parameter

    Returns:
        numpy.array: a 1-d vector of length n-1.

    Notes:
        For each probability density element two members of the input vectors y and t must be used. For example,
        if;
            x = [x_0, x_1, x_2] and
            t = [t_0, t_1, t_2]
        we can only obtain 2 probability density values;
            f_0 = f(x_1, t_1; x_0, t_0; a, b, sigma), the probability associated with evaluation point (x_1, t_1)
            f_1 = f(x_2, t_2,; x_1, t_1; a, b, sigma), the probability associated with the evaluation point (x_2, t_2)

        The returned result is the vector of length 2: [f_0, f_1]

    Doctests:
        >>> y = np.array([1.6,4.2,3,2.1])
        >>> t = np.array([1,2,3,4])
        >>> cirpdf(y, t, 1, 2, 0.2)
        array([  4.04336846e-23,   1.16777299e+00,   8.81371441e-01])
    """
    d = np.exp(-a * (t[1:] - t[:-1]))
    c = 2 * a / (sigma ** 2 * (1 - d))
    q = 2 * a * b / (sigma ** 2) - 1
    z = 2 * c * y[1:]  # <- transformed variable
    _lambda = 2 * c * y[:-1] * d  # <- non-centrality
    df = 2 * q + 2  # <- degrees of freedom
    if not (2 * a * b > sigma ** 2):  # <- check feller condition
        res = np.full(len(t) - 1, 1e-100)  # <- if fails return negligible probability density
    else:
        res = 2 * c * ncx2.pdf(z, df, _lambda)  # <- else return real probability density.
    return res
