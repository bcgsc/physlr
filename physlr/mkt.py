r"""
Python module to compute the Mann-Kendall test for trend in time series data.

This module contains a single function 'test' which implements the Mann-Kendall
test for a linear trend in a given time series. 

Introduction to the Mann-Kendall test
-------------------------------------

The Mann-Kendall test is used to determine whether or not there is a linear
monotonic trend in a given time series data. It is a non-parametric trend
closely related to the concept of Kendall's correlation coefficient [1]_. The
null hypothesis, :math:`H_0`, states that there is no monotonic trend, and this
is tested against one of three possible alternative hypotheses, :math:`H_a`:
(i) there is an upward monotonic trend, (ii) there is a downward monotonic
trend, or (iii) there is either an upward monotonic trend or a downward
monotonic trend. It is a robust test for trend detection used widely in
financial, climatological, hydrological, and environmental time series
analysis.

Assumptions underlying the Mann-Kendall test
--------------------------------------------

The Mann-Kendall test involves the following assumptions [2]_ regarding the 
given time series data:

    1. In the absence of a trend, the data are independently and identically
    distributed (iid).

    2. The measurements represent the true states of the observables at the
    times of measurements.

    3. The methods used for sample collection, instrumental measurements and
    data handling are unbiased.

Advantages of the Mann-Kendall test
-----------------------------------

The Mann-Kendall test provides the following advantages:

    1. It does not assume the data to be distributed according to any
    particular rule, i.e., e.g., it does not require that the data be normally
    distributed.

    2. It is not effected by missing data other than the fact the number of
    sample points are reduced and hence might effect the statistical
    significance adversely.

    3. It is not effected by irregular spacing of the time points of
    measurement.
    
    4. It is not effected by the length of the time series.

Limitations of the Mann-Kendall test
------------------------------------

The following limitations have to be kept in mind:
    1. The Mann-Kendall test is not suited for data with periodicities (i.e.,
    seasonal effects). In order for the test to be effective, it is recommended
    that all known periodic effects be removed from the data in a preprocessing
    step before computing the Mann-Kendall test.

    2. The Mann-Kendall test tends to give more negative results for shorter
    datasets, i.e., the longer the time series the more effective is the trend
    detection computation. 

Formulae
--------

The first step in the Mann-Kendall test for a time series :math:`x_1, x_2,
\dots, x_n` of length :math:`n` is to compute the indicator function
:math:`sgn(x_i - x_j)` such that:

    .. math::

        sgn(x_i - x_j) &= 
            \begin{cases}
                            1,  & x_i - x_j > 0\\
                            0,  & x_i - x_j = 0\\
                            -1, & x_i - x_j < 0
            \end{cases},

which tells us whether the difference between the measurements at time
:math:`i` and :math:`j` are positive, negative or zero.

Next, we compute the mean and variance of the above quantity. The mean
:math:`E[S]` is given by:

    .. math::

        E[S] = \sum_{i=1}^{n-1} \sum_{j=i+1}^{n} sgn(x_i - x_j),

and the variance :math:`VAR(S)` is given by:

    .. math::

        VAR(S) = \frac{1}{18} \Big( n(n-1)(2n+5) - \sum_{k=1}^p
        q_k(q_k-1)(2q_k+5) \Big),

where :math:`p` is the total number of tie groups in the data, and :math:`q_k`
is the number of data points contained in the :math:`k`-th tie group. For
example, if the time series measurements were {12, 56, 23, 12, 67, 45, 56, 56,
10}, we would have two tie groups for the measurements 12 and 56, i.e.
:math:`p=2`, and the number of data points in these tie groups would
:math:`q_1=2` for the tie group with {12}, and :math:`q_2=3` for the tie group
with {56}.

Using the mean :math:`E[S]` and the variance :math:`VAR(S)` we compute the
Mann-Kendall test statistic, using the following transformation, which ensures
that for large sample sizes, the test statistic :math:`Z_{MK}` is distributed
approximately normally:

    .. math::

        Z_{MK} &= 
            \begin{cases}
                            \frac{E[S] - 1} {\sqrt{VAR(S)}},  & E[S] > 0\\
                            0,  & E[S] = 0\\
                            \frac{E[S] + 1} {\sqrt{VAR(S)}},  & E[S] < 0\\
            \end{cases}.


Hypothesis testing
------------------

At a significance level :math:`\alpha` of the test, which is also the Type I
error rate, we compute whether or not to accept the alternative hypothesis 
:math:`H_a` for each variant of :math:`H_a` separately:

    :math:`H_a`: There exists an upward monotonic trend
        If :math:`Z_{MK} \geq Z_{1 - \alpha}` then accept :math:`H_a`, where the
        notation :math:`Z_{1 - \alpha}` denotes the :math:`100(1-\alpha)`-th
        percentile of the standard normal distribution.

    :math:`H_a`: There exists a downward monotonic trend
        If :math:`Z_{MK} \leq -Z_{1 - \alpha}` then accept :math:`H_a`.

    :math:`H_a`: There exists either an upward or a downward monotonic trend
        If :math:`|Z_{MK}| \geq Z_{1 - \alpha/2}` then accept :math:`H_a`,
        where the notation :math:`|\cdot|` is used to denote the absolute
        value function.

Updated formulae for implementation
-----------------------------------

One crucial notion involved in the Mann-Kendall test statistic is that of
whether the difference between two measurements is greater than, equal to, or
less than zero. This idea is in turn critically linked to the least count
(i.e., the minimum possible measurement value) of the time series measurements
:math:`x_i`. For example, let us consider the case when we measure :math:`x_i`
with a precision :math:`\varepsilon = 0.01`. In such a case, let us say for 
some reason, floating point errors in the entries of :math:`x_i` in the
memory, lead to a :math:`x_{11} - x_{27} = 0.000251 > 0`. However, to say that
this difference is actually greater than zero is meaningless! This is because
on the basis of the same measurement process we used on :math:`x`, we could
never ascertain such a small difference. This is why, in this implementation of
the Mann-Kendall test, we have included the least count error
:math:`\varepsilon` as a compulsory requirement for the test statistic
estimation.

This allows us to revise the above formulae fo rthe Mann-Kendall test as:

    .. math::

        sgn(x_i - x_j) &= 
            \begin{cases}
                            1,  & x_i - x_j > \varepsilon\\
                            0,  & |x_i - x_j| \leq \varepsilon\\
                            -1, & x_i - x_j < -\varepsilon
            \end{cases},

and:

    .. math::

        Z_{MK} &= 
            \begin{cases}
                            \frac{E[S] - 1} {\sqrt{VAR(S)}},  & E[S] >
                            \varepsilon\\
                            0,  & |E[S]| \leq \varepsilon\\
                            \frac{E[S] + 1} {\sqrt{VAR(S)}},  & E[S] <
                            -\varepsilon\\
            \end{cases}.

These revised formulae are the ones that are implemented in the :py:func:`test`
of this module.

Additional estimates
--------------------

In addition to the result of the Mann-Kendall test, which is in the form of a
string indicating whether or not to accept the alternative hypothesis, the
:py:func:`test` function also return a few additional estimates related to the
estimation of a monotonic trend in the time series.

Estimation of the simple linear regression parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The slope :math:`m` and intercept :math:`c` of a straight line fitted through
the time series data are estimated as follows:

    .. math::

        m = r_{x,t} \frac{\sigma_x}{\sigma_t},

where r_{x,t} is the Pearson's cross-correlation coefficient between
:math:`x` and :math:`t`.

    .. math::
        
        c = \mu_x - m \mu_t

where :math:`\mu` denotes the mean of the both variables respectively.

Estimation of :math:`p`-values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :py:func:`test` function also returns the :math:`p`-values for the given 
dataset under the various alternative hypotheses. Note that the estimation of
the :math:`p`-value is not essential to the computation of the test results as
formulated above. The :math:`p`-values need to estimated separately depending
on the type of alternative hypothesis used and the sign of :math:`E[S]`.
Denoting :math:`f(u)` as the probability density function of the standard
normal distribution, we can write down the :math:`p`-values as:

    :math:`H_a`: There exists an upward monotonic trend

    .. math::

        p_{Z_{MK}} &= 
            \begin{cases}
                \int_{Z_{MK}}^{\infty} f(u) \mathrm{d}u,& |E[S]|>\varepsilon\\
                0.5,  & |E[S]| \leq \varepsilon\\
            \end{cases}.

    :math:`H_a`: There exists a downward monotonic trend

    .. math::

        p_{Z_{MK}} &= 
            \begin{cases}
                \int^{Z_{MK}}_{-\infty} f(u) \mathrm{d}u,& |E[S]|>\varepsilon\\
                0.5,  & |E[S]| \leq \varepsilon\\
            \end{cases}.

    :math:`H_a`: There exists either an upward or a downward monotonic trend

    .. math::

        p_{Z_{MK}} &= 0.5 
            \begin{cases}
                \int_{Z_{MK}}^{\infty} f(u) \mathrm{d}u,& E[S]>\varepsilon\\
                1,  & |E[S]| \leq \varepsilon\\
                \int^{Z_{MK}}_{-\infty} f(u) \mathrm{d}u,& E[S]<-\varepsilon\\
            \end{cases}.


References
----------
.. [1]  | Pohlert, T.
        | "Non-Parametric Trend Tests and Change-Point Detection".
        | R-package `trend`. Accessed on: 17 April, 2017.
        | https://cran.r-project.org/web/packages/trend/vignettes/trend.pdf

.. [2]  | "Mann-Kendall Test For Monotonic Trend".
        | Visual Simple Plan. Accessed on: 17 April, 2017.
        | http://vsp.pnnl.gov/help/Vsample/Design_Trend_Mann_Kendall.htm

"""
# Created: Mon Apr 17, 2017  01:18PM
# Last modified: Mon Apr 17, 2017  09:24PM
# Copyright: Bedartha Goswami <goswami@uni-potsdam.de>


import numpy as np
from scipy.special import ndtri, ndtr
import sys


def test(t, x, eps=None, alpha=None, Ha=None):
    """
    Runs the Mann-Kendall test for trend in time series data.

    Parameters
    ----------
    t : 1D numpy.ndarray
        array of the time points of measurements
    x : 1D numpy.ndarray
        array containing the measurements corresponding to entries of 't'
    eps : scalar, float, greater than zero
        least count error of measurements which help determine ties in the data
    alpha : scalar, float, greater than zero
        significance level of the statistical test (Type I error)
    Ha : string, options include 'up', 'down', 'upordown'
        type of test: one-sided ('up' or 'down') or two-sided ('updown')

    Returns
    -------
    MK : string
        result of the statistical test indicating whether or not to accept hte
        alternative hypothesis 'Ha'
    m : scalar, float
        slope of the linear fit to the data
    c : scalar, float
        intercept of the linear fit to the data
    p : scalar, float, greater than zero
        p-value of the obtained Z-score statistic for the Mann-Kendall test

    Raises
    ------
    AssertionError : error
                    least count error of measurements 'eps' is not given
    AssertionError : error
                    significance level of test 'alpha' is not given
    AssertionError : error
                    alternative hypothesis 'Ha' is not given

    """
    # assert a least count for the measurements x
    assert eps, "Please provide least count error for measurements 'x'"
    assert alpha, "Please provide significance level 'alpha' for the test"
    assert Ha, "Please provide the alternative hypothesis 'Ha'"

    # estimate sign of all possible (n(n-1)) / 2 differences
    n = len(t)
    sgn = np.zeros((n, n), dtype="int")
    for i in range(n):
        tmp = x - x[i]
        tmp[np.where(np.fabs(tmp) <= eps)] = 0.
        sgn[i] = np.sign(tmp)

    # estimate mean of the sign of all possible differences
    S = sgn[np.triu_indices(n, k=1)].sum()

    # estimate variance of the sign of all possible differences
    # 1. Determine no. of tie groups 'p' and no. of ties in each group 'q'
    np.fill_diagonal(sgn, eps * 1E6)
    i, j = np.where(sgn == 0.)
    ties = np.unique(x[i])
    p = len(ties)
    q = np.zeros(len(ties), dtype="int")
    for k in range(p):
        idx =  np.where(np.fabs(x - ties[k]) < eps)[0]
        q[k] = len(idx)
    # 2. Determine the two terms in the variance calculation
    term1 = n * (n - 1) * (2 * n + 5)
    term2 = (q * (q - 1) * (2 * q + 5)).sum()
    # 3. estimate variance
    varS = float(term1 - term2) / 18.

    # Compute the Z-score based on above estimated mean and variance
    if S > eps:
        Zmk = (S - 1) / np.sqrt(varS)
    elif np.fabs(S) <= eps:
        Zmk = 0.
    elif S < -eps:
        Zmk = (S + 1) / np.sqrt(varS)

    # compute test based on given 'alpha' and alternative hypothesis
    # note: for all the following cases, the null hypothesis Ho is:
    # Ho := there is no monotonic trend
    # 
    # Ha := There is an upward monotonic trend
    if Ha == "up":
        Z_ = ndtri(1. - alpha)
        if Zmk >= Z_:
            MK = "accept Ha := upward trend"
        else:
            MK = "reject Ha := upward trend"
    # Ha := There is a downward monotonic trend
    elif Ha == "down":
        Z_ = ndtri(1. - alpha)
        if Zmk <= -Z_:
            MK = "accept Ha := downward trend"
        else:
            MK = "reject Ha := downward trend"
    # Ha := There is an upward OR downward monotonic trend
    elif Ha == "upordown":
        Z_ = ndtri(1. - alpha / 2.)
        if np.fabs(Zmk) >= Z_:
            MK = "accept Ha := upward OR downward trend"
        else:
            MK = "reject Ha := upward OR downward trend"

    # ----------
    # AS A BONUS
    # ----------
    # estimate the slope and intercept of the line
    m = np.corrcoef(t, x)[0, 1] * (np.std(x) / np.std(t))
    c = np.mean(x) - m * np.mean(t)

    # ----------
    # AS A BONUS
    # ----------
    # estimate the p-value for the obtained Z-score Zmk
    if S > eps:
        if Ha == "up":
            p = 1. - ndtr(Zmk)
        elif Ha == "down":
            p = ndtr(Zmk)
        elif Ha == "upordown":
            p = 0.5 * (1. - ndtr(Zmk))
    elif np.fabs(S) <= eps:
        p = 0.5
    elif S < -eps:
        if Ha == "up":
            p = 1. - ndtr(Zmk)
        elif Ha == "down":
            p = ndtr(Zmk)
        elif Ha == "upordown":
            p = 0.5 * (ndtr(Zmk))

    return MK, m, c, p
