from sage.all import *
import Crypto.Util.number as number 

def sign(x):
    """
    Returns the sign of input variable [x].
    """
    if x == 0:
        return 0
    elif x < 0:
        return -1
    else:
        return 1

def long_to_bytes(x):
    return number.long_to_bytes(int(x))

def bytes_to_long(x):
    return Integer(number.bytes_to_long(x))

def byte_length(x):
    return ceil(log(x, 256))

def find_roots_univariate(x, polynomial):
    """
    Yields all roots of a univariate polynomial [polynomial] in variable [x].
    """
    if polynomial.is_constant():
        return

    for root in polynomial.roots(multiplicities=False):
        if root != 0:
            yield {x: int(root)}

def find_common_root_bivariate(polynomial_ring, poly1, poly2):
    """
    Yields a common root of two bivariate polynomial [poly1, poly2] in polynomial ring [polynomial_ring] by Groebner basis method.
    """
    assert polynomial_ring.ngens() == 2

    I = Sequence([poly1, poly2], polynomial_ring.change_ring(QQ, order="lex")).ideal()
    G = I.groebner_basis()

    roots = {}
    for polynomial in G:
        if len(polynomial.variables()) == 1:
            x = polynomial.variables()[0]
            for root in find_roots_univariate(x, polynomial.univariate_polynomial()):
                roots |= root

    if roots == {}:
        return
    if len(roots) == 2:
        yield roots
        return

    for polynomial in G:
        if polynomial.variables() != [] or polynomial.variables() != [x]:
            subpoly = polynomial.subs(roots)
            y = subpoly.variables()[0]
            for root in find_roots_univariate(y, subpoly.univariate_polynomial()):
                roots |= root
    if len(roots) == 2:
        yield roots
        return