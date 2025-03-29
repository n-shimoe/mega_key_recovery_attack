import Crypto.Util.number as number
from Crypto.Cipher import AES
from time import process_time
from shared.util import *
from shared.constant import *
from shared.mega_sim import *

# Parameters of PoC
iteration_number = 100
b = 4
t = b-2
print("b:", b)
print("t:", t)
assert t <= b-2 and t <= 6 and t >= 0

# Variables to store measured values
success_number = 0
runtime_list = []
lattice_dimension_list = []
gamma1_list = []
gamma2_list = []
quotient_number_list = []

def recover_q_by_solving_approximate_divisor(a, N, beta, X):
    """
    Solves the approximate divisor problem:
        [a] + x = 0 (mod q) (|x| <= [X]),
    where
    - q is an unknown divisor of [N],
    - q >= [N]^[beta],
    - q = [a] + y (y is a solution of the above equation),
    and return q if the recovery of q is successful or None otherwise.
    """
    global lattice_dimension_list

    gamma = log(X, N).n()
    m_ad = ceil((1 - beta)/(beta^2 - gamma)) - 1
    t_ad = ceil(beta * m_ad)
    Ring.<x> = PolynomialRing(ZZ, order="lex")

    polylist = []
    for i in range(m_ad+1):
        f = (a + x)^i * N^(max(t_ad - i, 0))
        polylist.append(f(X*x))
    dim = len(polylist)
    print("Lattice dimension:", dim)
    lattice_dimension_list.append(dim)
    
    coeffs = []
    for poly in polylist:
        row = []
        for i in range(m_ad+1):
            row.append(poly.monomial_coefficient(x^i))
        coeffs.append(row)

    B = Matrix(ZZ, coeffs)
    reduced_B = B.LLL()

    polynomial = 0*x
    for j in range(dim):
        polynomial += reduced_B[0][j] * x^j / X^j

    roots = find_roots_univariate(x, polynomial)
    for root in roots:
        y = int(root[x])
        estimated_q = gcd(a + y, N)
        if estimated_q != N and estimated_q != 1:
            return estimated_q
    return None

def recover_q_by_solving_ad_with_help_of_approximation_of_q2(a1, a2, N, beta, X1, X2):
    """
    Recovers an unknown divisor q of N under the following condition.
    - a1 + e1 = q (|e1|<=X1)
    - a2 + e2 = q^2 (|e2|<=X2)
    If its recovery fails, return None.
    """
    a1 += X1
    X1 *= 2

    alpha = ceil((a2 - a1^2)/(2*a1))
    refined_a1 = a1 + alpha
    refined_X1 = ceil(max(X2, X1^2) / a1)

    return recover_q_by_solving_approximate_divisor(refined_a1, N, beta, refined_X1)

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

    start_time = process_time()
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
        q = recover_q_by_solving_ad_with_help_of_approximation_of_q2(q_approx, q2_approx, N, 0.5, q_error_bound, q2_error_bound)
        if q != None:
            success_number += 1
            print("Prime q is recovered")
            break
    
    end_time = process_time()
    runtime = end_time - start_time
    runtime_list.append(runtime)
    print("Runtime (s):", runtime)

# Main
for iteration in range(iteration_number):
    print(f"Iteration {iteration+1}")
    attack_simluated_MEGA(b, t)
    print()

print("Success rate (%):", 100*success_number/iteration_number)
print("Average of γ1:", sum(gamma1_list).n()/len(gamma1_list))
print("Average of γ2:", sum(gamma2_list).n()/len(gamma2_list))
print("Average of runtime (s):", sum(runtime_list)/len(runtime_list))