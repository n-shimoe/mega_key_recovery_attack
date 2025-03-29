from Crypto.Cipher import AES
from sage.matrix.constructor import Matrix
import logging
from time import process_time
from shared.mega_sim import *
from shared.util import *

# Parameters of PoC
iteration_number = 10
b = 4
t = b-2
print("b:", b)
print("t:", t)
assert t <= b-2 and t <= 6 and t >= 0

# Variables to store measured values
success_number = 0
runtime_list = []
lattice_dimension_list = []
ideal_dimension_list = []
gamma1_list = []
gamma2_list = []
quotient_number_list = []

def solve_simultaneous_modular_equations(a1, a2, r1, r2, N, beta, X1, X2):
    """
    Solves the followng simultaneous modular equations:
        [a1] + x1 = 0 (mod q^[r1]) 
        [a2] + x2 = 0 (mod q^[r2]),
    where q is an unknown divisor of [N], q >= N^[beta], |x1| <= [X1], and |x2| <= [X2].
    """
    global lattice_dimension_list
    global ideal_dimension_list

    t = 0
    g1 = log(X1, N).n()
    g2 = log(X2, N).n()
    # Determine the parameter t for the algorithm
    for m in range(20):
        w = 0
        s_N = 0
        s_X1 = 0
        s_X2 = 0
        for i1 in range(ceil(m/g1)):
            for i2 in range(ceil(m/g2)):
                if g1*i1 + g2*i2 <= beta*m:
                    w += 1
                    s_X1 += i1
                    s_X2 += i2
                if r1*i1 + r2*i2 <= m:
                    s_N += ceil(m - r1*i1 - r2*i2)
        if g1*s_X1 + g2*s_X2 + s_N < beta*m*(w-1):
            t = m
            break
    if t == 0:
        print("Too large lattice dimension to solve")
        return (0, 0)

    R.<x1, x2> = PolynomialRing(ZZ)
    polylist = []
    index = []
    for i1 in range(ceil(t/g1)):
        for i2 in range(ceil(t/g2)):
            if g1*i1 + g2*i2 <= beta*t:
                f = (a1 + x1)^i1 * (a2 + x2)^i2 * N^(max(t - r1*i1 - r2*i2, 0))
                polylist.append(f(X1*x1, X2*x2))
                index.append((i1, i2))

    lattice_dim = len(polylist)
    print(f"Lattice dimension:", lattice_dim)
    lattice_dimension_list.append(lattice_dim)

    coeffs = []
    for poly in polylist:
        row = []
        for i1, i2 in index:
            row.append(poly.coefficient({x1: i1, x2: i2}))
        coeffs.append(row)

    B = Matrix(ZZ, coeffs)
    reduced_B = B.LLL()

    print("Howgrave-Graham's lemma is satisfied?:", ceil(reduced_B[1].norm())  * ceil(sqrt(lattice_dim)) < N^(ceil(beta*t)))

    short_polys = []
    for i in range(lattice_dim):
        h = 0*x1
        for j, ind in enumerate(index):
            i1, i2 = ind
            h += reduced_B[i][j] / (X1^i1 * X2^i2) * x1^i1 * x2^i2
        short_polys.append(h)
    
    ideal_dim = Sequence(short_polys[:2], R.change_ring(QQ, order="lex")).ideal().dimension()
    ideal_dimension_list.append(ideal_dim)
    print("Ideal dimension:", ideal_dim)

    roots = find_common_root_bivariate(R, short_polys[0], short_polys[1])
    for root in roots:
        yield (root[x1], root[x2])

def recover_q_by_solving_simultaneous_modular_equations(a1, a2, N, beta, X1, X2):
    """
    Tries to recover an unknown divisor q of [N] as shown below:
        [a1] + x1 = q
        [a2] + x2 = q^2,
    where q >= N^[beta], |x1| <= [X1], and |x2| <= [X2].
    """
    solutions = solve_simultaneous_modular_equations(a1, a2, 1, 2, N, 0.5, X1, X2)
    for sol in solutions:
        y1, y2 = sol
        estimated_q = gcd(q_approx + y1, N)
        if estimated_q != N and estimated_q != 1:
            return q
    return None


def attack_simluated_MEGA(b, t):
    global success_number
    global gamma1_list
    global gamma2_list

    # Skip the parts of solving HNP-SUM and recovering the MSB of q 
    k, approx_hnp_sum, N, q_M, q_approx, q_error_bound = get_oracle_output(b, t)
    q2_error_bound = ceil((1 << HNP_SUM_ERROR_SIZE)/k)
    gamma1 = log(q_error_bound, N).n()
    gamma2 = log(q2_error_bound, N).n()
    gamma1_list.append(gamma1)
    gamma2_list.append(gamma2)
    print("γ1:", gamma1)
    print("γ2:", gamma2)

    # Enumerate the candidates of quotient of k*q^2
    Quo = (k * q_M^2 * (1 << (2*(PRIME_SIZE - AES_BLOCK_SIZE*b + 16)))) // N
    R = (k * q_M^2 * (1 << (2*(PRIME_SIZE - AES_BLOCK_SIZE*b + 16)))) % N
    Ep = q_M * (1 << (2*(PRIME_SIZE - AES_BLOCK_SIZE*b + 16)+2))
    Quo_dd = (R + k*Ep)//N + 1
    quotient_list = [Quo + l for l in range(Quo_dd)]
    print("The number of quotient candidates:", len(quotient_list))
    quotient_number_list.append(len(quotient_list))

    for Q in quotient_list:
        # Create the approximated value of q^2 from a HNP sample
        q2_approx = ceil((Q*N + approx_hnp_sum)/k)

        # Recover q
        q = recover_q_by_solving_simultaneous_modular_equations(q_approx, q2_approx, N, 0.5, q_error_bound, q2_error_bound)
        if q != None:
            success_number += 1
            print("Prime q is recovered")
            break

# Main
for iteration in range(iteration_number):
    print(f"Iteration {iteration+1}")
    attack_simluated_MEGA(b, t)
    print()    
    
print("Success rate (%):", 100*success_number/iteration_number)
print("Average of γ1:", sum(gamma1_list).n()/len(gamma1_list))
print("Average of γ2:", sum(gamma2_list).n()/len(gamma2_list))
print("Average of lattice dimension: ", sum(lattice_dimension_list)/len(lattice_dimension_list))
print("Average of ideal dimension: ", sum(ideal_dimension_list)/len(ideal_dimension_list))