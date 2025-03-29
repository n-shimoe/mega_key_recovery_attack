from Crypto.Cipher import AES
from sage.all import *
from shared.constant import *
from shared.util import *

def create_RSA_instance():
    """
    Returns CRT-RSA instance [e, N, p, q, u, d].
    - Public key (e, N): N=p*q.
    - Secret key (p, q, d, u): p, q are distinct primes. d = e^{-1} (mod (p-1)(q-1)). u = q^{-1} (mod p).
    Assume that the byte length of d, u is 256, 128, respectively.
    """
    is_desired_instance = False
    while not is_desired_instance:
        e = PUBLIC_EXPONENT
        p = next_prime(ZZ.random_element(1 << (PRIME_SIZE-1), 1 << PRIME_SIZE))
        while gcd(e, p-1) != 1:
            p = next_prime(ZZ.random_element(1 << (PRIME_SIZE-1), 1 << PRIME_SIZE))
        q = next_prime(ZZ.random_element(1 << (PRIME_SIZE-1), 1 << PRIME_SIZE))
        while gcd(e, q-1) != 1 or p == q:
            q = next_prime(ZZ.random_element(1 << (PRIME_SIZE-1), 1 << PRIME_SIZE))
        d = Integer(pow(e, -1, (p-1)*(q-1)))
        u = Integer(pow(q, -1, p))
        N = p*q
        if byte_length(d) == 256 and byte_length(u) == 128:
            is_desired_instance = True

    return e, N, p, q, u, d

def decrypt_rsa(c, p, q, u, d):
    """
    Decrypts [c] under the given CRT-RSA secret key ([p], [q], [d], [u]),
    and returns the result of decryption.
    """
    dp = d % (p-1)
    dq = d % (q-1)
    mp = Integer(pow(c, dp, p))
    mq = Integer(pow(c, dq, q))
    return ( ((((mp - mq) % p) * u) % p) * q + mq ) % (p*q)

def encode_rsa_secret_key(p, q, d, u):
    encoded_p = long_to_bytes(p)
    encoded_q = long_to_bytes(q)
    encoded_u = long_to_bytes(u)
    encoded_d = long_to_bytes(d)
    
    length_p = long_to_bytes(p.nbits())
    length_q = long_to_bytes(q.nbits())
    length_u = long_to_bytes(u.nbits())
    if len(length_u) < 2:
        length_u = b'\x00' + length_u
    length_d = long_to_bytes(d.nbits())
    if len(length_d) < 2:
        length_d = b'\x00' + length_d
    concatenated_sk = length_q + encoded_q + length_p + encoded_p + length_d + encoded_d + length_u + encoded_u
    
    # add padding P to concatenated_sk
    while len(concatenated_sk) % 16 != 0:
        concatenated_sk += b'\x00' # Note: random padding in the actual MEGA
    
    return concatenated_sk

def decode_rsa_secret_key(encoded_rsa_secret_key):
    """
    Decodes the given encoded RSA secret key [encode_rsa_secret_key].
    """
    rsa_secret_key = []
    index = 0
    for _ in range(4):
        length = ceil(bytes_to_long(encoded_rsa_secret_key[index:index+2])/8)
        index = index+2
        value = bytes_to_long(encoded_rsa_secret_key[index:index+length])
        rsa_secret_key.append(value)
        index = index + length
    return rsa_secret_key

def get_block(block_index, aes_byte_string):
    """
    Obtains a block of AES plaintext or ciphertext [aes_byte_string] corresponding to the designated index [block_index].
    Note that the indexing is zero-based.
    """
    if block_index >= 0 and block_index < len(aes_byte_string)//(AES_BLOCK_SIZE//8):
        start = block_index * (AES_BLOCK_SIZE//8)
        end = (block_index + 1) * (AES_BLOCK_SIZE//8)
        desired_block = aes_byte_string[start:end]
        return desired_block
    else:
        print("Get_block process is failed! Returned value is empty string.")
        return b""

def overwrite_block(block_index, aes_byte_string, new_block):
    """
    Overwrites a block of AES plaintext or ciphertext [aes_byte_string] corresponding to the designated index [block_index]
    with a block [new_block].
    Note that the indexing is zero-based.
    """
    if block_index >= 0 and block_index < len(aes_byte_string)//(AES_BLOCK_SIZE//8):
        if len(new_block)*8 != AES_BLOCK_SIZE:
            print("Overwrite process is failed because of invalid block size! Returned value is the original aes_byte_string.")
            return aes_byte_string
        else:
            start = block_index * (AES_BLOCK_SIZE//8)
            end = (block_index + 1) * (AES_BLOCK_SIZE//8)
            new_aes_byte_string = aes_byte_string[:start] + new_block + aes_byte_string[end:]
            return new_aes_byte_string
    else:
        print("Overwrite process is failed because of invalid block index! Returned value is the original aes_byte_string.")
        return aes_byte_string

def optimal_overwriting_for_rsa_secret_key(encrypted_rsa_secret_key, t):
    block_index = [17, 0, 1, 2, 3]
    blocks = []
    for i in block_index:
        blocks.append(get_block(i, encrypted_rsa_secret_key))

    overwritten_secret_keys = []
    for block in blocks:
        overwritten_secret_key = encrypted_rsa_secret_key
        overwritten_secret_key = overwrite_block(39-t, overwritten_secret_key, block)
        overwritten_secret_keys.append(overwritten_secret_key)
    
    return overwritten_secret_keys

def get_oracle_output(b, t):
    """
    Outputs (k, a, N, q_M, a_q, X) for the number [b]+1 of login attempts
    and the overwriting position 39-[t].
    Note that the indexing is zero-based.

    RSA modulus: N = p*q
    MSBs of RSA prime q: q_M
    HNP-SUM sample: delta*x = a + e (mod N), (|e| <= 2^HNP_SUM_ERROR_SIZE, delta > 0, k=delta*2^(AES_BLOCK_SIZE*t+64))
    approximate divisor problem: a_q + y = 0 (mod q), (|y| <= X, a_q = q_M * 2^(PRIME_SIZE-AES_BLOCK_SIZE*b+16))
    """

    # RSA
    e, N, p, q, u, d = create_RSA_instance()

    # AES-ECB
    master_key = long_to_bytes(ZZ.random_element(1 << (AES_BLOCK_SIZE-1), 1 << AES_BLOCK_SIZE))
    aes_cipher = AES.new(master_key, AES.MODE_ECB)

    # Key overwriting
    encoded_sk = encode_rsa_secret_key(p, q, d, u)
    encrypted_encoded_sk = aes_cipher.encrypt(encoded_sk)
    overwritten_secret_keys = optimal_overwriting_for_rsa_secret_key(encrypted_encoded_sk, t)

    # Encrypt a message with RSA
    message = (1 << PRIME_SIZE)
    c = pow(message, e, N)

    # Create HNP-SUM instance
    approximation_samples = []
    error_list = [] # For test
    for overwritten_sk in overwritten_secret_keys:
        decrypted_encoded_sk = aes_cipher.decrypt(overwritten_sk)
        decoded_sk = decode_rsa_secret_key(decrypted_encoded_sk)
        q, p, d, u_dash = decoded_sk
        decryption_result = decrypt_rsa(c, p, q, u_dash, d)
        assert message != decryption_result
        a = ((decryption_result >> HNP_SUM_ERROR_SIZE) << (HNP_SUM_ERROR_SIZE)) # approximated value
        err = (decryption_result % (1 << HNP_SUM_ERROR_SIZE))
        approximation_samples.append(a)
        error_list.append(err)

    # Get MSB of q, and create the approximated valeu of q
    q_M = q >> (PRIME_SIZE - AES_BLOCK_SIZE*b + 16)
    q_approx = q_M << (PRIME_SIZE - AES_BLOCK_SIZE*b + 16)
    q_error_bound = 1 << (PRIME_SIZE - AES_BLOCK_SIZE*b + 16)

    # Get answer of HNP-SUM
    block1 = bytes_to_long(get_block(17, encoded_sk))
    block2 = bytes_to_long(get_block(0, encoded_sk))
    delta = block2 - block1
    assert abs(delta) <= (1 << AES_BLOCK_SIZE)
    approx = approximation_samples[1] - approximation_samples[0]
    error = error_list[1] - error_list[0]
    k = delta << (AES_BLOCK_SIZE*t + 64)

    # Make delta positive
    delta_sign = sign(delta)
    k *= delta_sign
    approx = (delta_sign * approx) % N
    error *= delta_sign

    assert abs((k*q*q) % N - approx) <  (1 << HNP_SUM_ERROR_SIZE) # This HNP-SUM instance follows Heninger-Ryan's model of MEGA

    return k, approx, N, q_M, q_approx, q_error_bound