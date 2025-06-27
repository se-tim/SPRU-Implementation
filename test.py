from sage.all import log, randint
from spru_package.ckks_x import CKKS_x, Poly

# Configuration

if input("Enter your own parameters (y / n)? ") == "y":
    N = 2 ** int(input("Ring degree: N = 2^"))
    C = 2 ** int(input("Number of coefficients: C = 2^"))
    L = int(input("Maximal level during bootstrapping: L = "))
    q = 2 ** int(input("Base modulus: q = 2^"))
    p = 2 ** int(input("Scaling factor outside of bootstrapping: p = 2^"))
    delta = 2 ** int(input("Scaling factor during bootstrapping: delta = 2^"))
    h = 2 ** int(input("Hamming weight of secret key: h = 2^"))
    g = int(input("Grouping parameter (log_2 of radix) for SlotToCoeff: g = "))
else:
    N = 2**15
    C = 2**2
    L = 9
    q = 2**28
    p = 2**22
    delta = 2**35
    h = 2**6
    g = 1
    print("We chose the following parameters:")
    print(f"Ring degree: N = 2^{log(N, 2)}")
    print(f"Number of coefficients: C = 2^{log(C, 2)}")
    print(f"Maximal level during bootstrapping: L = {L}")
    print(f"Base modulus: q = 2^{log(q, 2)}")
    print(f"Scaling factor outside of bootstrapping: p = 2^{log(p, 2)}")
    print(f"Scaling factor during bootstrapping: delta = 2^{log(delta, 2)}")
    print(f"Hamming weight of secret key: h = 2^{log(h, 2)}")
    print(f"Grouping parameter (log_2 of radix) for SlotToCoeff: g = {g}")
print()

CKKS_x.config(N, L, q, p, delta, print_messages=True)
print()

CKKS_x.key_gen(h=h, print_messages=True)
print()

CKKS_x.config_SPRU(CKKS_x.sk, C, g, print_messages=True)
print()

print("Estimated security in bits:")
print(CKKS_x.get_security(check_primal_hybrid=False))
print()

# Testing SPRU bootstrapping

print("A randomly generated ciphertext:")
ct = CKKS_x.get_random_ciphertext(CKKS_x.sk)
print(ct)
print()

print("It decrypts to:")
print(CKKS_x.decode(ct.dec_to_poly(CKKS_x.sk), C // 2))
print()

print("Projecting the ciphertext down to the lowest level:")
ct_low = ct % q
print(ct_low)
print()

print("Bootstrapping it:")
ct_boot = ct_low.SPRU()
print(ct_boot)
print()

print("This decrypts to:")
print(CKKS_x.decode(ct_boot.dec_to_poly(CKKS_x.sk), C // 2))
print()

print("Precision in bits:")
print(ct_low.get_precision(ct_boot, CKKS_x.sk, C // 2))
