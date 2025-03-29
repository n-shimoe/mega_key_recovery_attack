import Crypto.Util.number as number
from Crypto.Cipher import AES
from time import process_time
from shared.util import *
from shared.constant import *
from shared.mega_sim import *

# Parameters of PoC
iteration_number = 100
b = 5
t = 0
print("b:", b)
print("t:", t)
assert b >= 5 and t <= 6 and t >= 0

# Variables to store measured values
success_number = 0
runtime_list = []
lattice_dimension_list = []
gamma1_list = []

def recover_q_by_Coppersmith(a, N, X):
    """
    Recovers q satisfying the following condition by Coppersmith's method in Heninger-Ryan's attack.
    - q is an unknown divisor of RSA modulus [N]
    - |q - [a]| <= [X]
    If its recovery fails, return None.
    """
    global lattice_dim

    R.<x> = PolynomialRing(ZZ)

    m = 4
    k = 2

    polylist = []
    for i in range(m+1):
        f = N^(min(k, m-i)) * (a+x)^i
        polylist.append(f(X*x))
    lattice_dim = m+1
    print(f"Lattice dimension:", lattice_dim)
    lattice_dimension_list.append(lattice_dim)
        
    coeffs = []
    for poly in polylist:
        row = []
        for i in range(m+1):
            row.append(poly.monomial_coefficient(x^i))
        coeffs.append(row)

    B = Matrix(ZZ, coeffs)
    reduced_B = B.LLL()

    polynomial = 0*x
    for j in range(m+1):
        polynomial += reduced_B[0][j] * x^j / X^j

    roots = find_roots_univariate(x, polynomial)
    for root in roots:
        y = int(root[x])
        estimated_q = gcd(a + y, N)
        if estimated_q != N and estimated_q != 1:
            return estimated_q
    return None

def attack_simluated_MEGA(b, t):
    global success_number
    global gamma1_list

    # Skip the parts of solving HNP-SUM and recovering the MSB of q
    _, _, N, _, q_approx, q_error_bound = get_oracle_output(b, t)
    gamma1 = log(q_error_bound, N).n()
    gamma1_list.append(gamma1)
    print("γ1:", gamma1)

    # Recover q by Coppersmith
    start_time = process_time()
    q = recover_q_by_Coppersmith(q_approx, N, q_error_bound)
    if q != None:
        assert gcd(q, N) != 1 and gcd(q, N) != N
        success_number += 1
        print("Prime q is recovered")
    end_time = process_time()
    runtime = end_time - start_time
    runtime_list.append(runtime)
    print("Runtime (s):", runtime)

# Main
for iteration in range(iteration_number):
    print(f"Iteration {iteration+1}")
    attack_simluated_MEGA(b, t)
    print()

print("Success rate(%):", 100*success_number/iteration_number)
print("Average of γ1:", sum(gamma1_list).n()/len(gamma1_list))
print("Average of lattice dimension: ", sum(lattice_dimension_list)/len(lattice_dimension_list))
print("Average of runtime:", sum(runtime_list)/len(runtime_list))