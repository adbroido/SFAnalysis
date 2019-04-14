import mpmath as mp
import numpy as np

""" Contains functions to compute the normalization constants for both power-law
and power-law with expoenential cutoff. Both functions are meant to be called
directly.

Note that plwcconstant returns a constant to divide by, while plconst returns
a constant to multiply.

"""



def plwcconst(alpha, lam, xmin):
    """ Computes the normalization constant on the discrete power law;
    i.e., computes C so that

        1/C * sum from xmin to infinity of ( x^(-alpha) * e^(-lam x)] ) = 1

    The formula below is obtained by noting that the sum above is equal to

        sum from xmin to infinity of ( x^(-alpha) * e^(-lam x)] ) =
        e^(-xmin * lam) * HurwitzLerchPhi(e^(-lam), alpha, xmin)

    where HurwitzLerchPhi is the Lerch Phi function:

        Phi(z,s,a) = sum from 0 to infinity of ( z^k / (a+k)^s )

    (as defined in sympy docs). If one is disinclined to note this fact
    by deriving it, one is encouraged to look at it in Mathematica and
    blindly trust the result.

    Note: The standard notation for the Lerch Phi function above uses "a" as
          a parameter. This does not correspond to the alpha we pass to the
          function.

    Inputs:
        alpha                  float, exponent on x, must be > -1
        lam                    float, exponential cutoff, must be > 0
        xmin                   int, starting point for sum, must be >= 1

    Outputs:
        C                      float, normalization constant

    """
    mp.mp.dps = 40
    result = mp.exp(-xmin * lam) * mp.lerchphi(mp.exp(-lam), alpha, xmin)
    C = float(result)
    return C

def plconst(alpha, xmin):
    """ Computes the normalization constant on the discrete power law;
    i.e., computes C so that

        C * sum from xmin to infinity of x^(-alpha) = 1

    The formula below is obtained by noting that the sum above is equal to

        sum from xmin to infinity of x^(-aplha) = HurwitzZeta(alpha, xmin)

    where HurwitzZeta is defined by

        zeta(s,a) = sum from 0 to infinity of 1 / (a+k)^s

    (as defined in sympy docs). If one is disinclined to note this fact
    by deriving it, one is encouraged to look at it in Mathematica and
    blindly trust the result.

    Note: The standard notation for the zeta function above uses "a" as
          a parameter. This does not correspond to the alpha we pass to the
          function. (Our alpha is s in the notation above).


    Inputs:
        alpha                  array, shape=(1,) exponent on x, must be > 1
                               must be passed as array for op.minimize()
        xmin                   int, starting point for sum, must be >= 1

    Outputs:
        C                      float, normalization constant

    """
    total = mp.zeta(np.asscalar(alpha),1) # op.minimize passes array
    lowertail = np.sum(np.asarray(range(1,xmin)**(-alpha)))
    result = total-lowertail
    C = 1./(result)
    return float(C)



# example usage:
if __name__ == '__main__':
    alpha = np.array([1.2])
    lam = 2.0
    xmin = 3

    const1 = plwcconst(alpha, lam, xmin)
    const2 = plconst(alpha, xmin)
