from sage.all import log, randint
import numpy as np
from .ckks_in_sagemath.ckks_package.ckks import CKKS
from .ckks_in_sagemath.ckks_package.fast_dft import get_grouped_F
from .ckks_in_sagemath.ckks_package.poly import Poly
from .ext_bit_rev import ext_bit_rev_vector


class CKKS_x(CKKS):
    """
    Class for CKKS scheme with SPRU bootstrapping.
    """

    @classmethod
    def config(cls, N, L_boot, q0, p, delta, print_messages=False):
        """
        Configure basic CKKS parameters.

        Args:
            N (int):
                Ring dimension (must be a power of two).
            L_boot (int):
                Maximal level during bootstrapping.
            q0 (int):
                Smallest modulus.
            p (int):
                Scaling factor outside of bootstrapping.
            delta (int):
                Scaling factor during bootstrapping.
            print_messages (bool, optional):
                Whether to print information messages. Defaults to False.
        """
        super().config(
            N, N // 2, L_boot, q0, p, delta, print_messages=print_messages
        )

    @classmethod
    def key_gen(cls, h=None, print_messages=False):
        """
        Generate the secret key, public key and evaluation key for the scheme.

        Args:
            h (int, optional):
                Hamming weight of the secret key. Must be a power of two.
                Defaults to 2 ** min(6, (log(N, 2) // 2)).
            print_messages (bool, optional):
                Whether to print information messages. Defaults to False.
        """

        if h is not None and h & (h - 1) != 0:
            raise ValueError("h must be a power of two.")
        if h is None:
            h = 2 ** min(6, int(log(cls.N, 2) // 2))

        # Generate the SPRU secret key

        B = cls.N // h
        sk_coeffs = np.zeros(cls.N, dtype=int)
        sk_coeffs[0] = 1
        for b in range(1, h):
            j = randint(0, B - 1)
            sk_coeffs[b * B + j] = 1
        sk = Poly(sk_coeffs, cls.N)

        q = cls.moduli_boot[-2]  # Largest modulus for evaluation key
        P = cls.moduli_boot[-2]  # Extra factor added to evaluation key modulus
        super().key_gen(sk=sk, q=q, P=P, print_messages=print_messages)

    @classmethod
    def get_security(cls, check_primal_hybrid=False):
        """
        Use the LWE estimator to compute the security level.

        Args:
            check_primal_hybrid (bool, optional):
                Whether to check the primal hybrid security level. Defaults to
                False.

        Returns:
            float:
                Logarithm (base 2) of estimated number of operations required
                to break the scheme.
        """
        return super().get_security(
            N=cls.N - cls.N // cls.h,
            q=cls.moduli_boot[-2] ** 2,
            sk_plus=cls.h - 1,
            sk_minus=0,
            check_primal_hybrid=check_primal_hybrid,
        )

    @classmethod
    def config_SPRU(cls, sk, C, g, print_messages=False):
        """
        Perform precomputations needed for SPRU bootstrapping.

        Args:
            sk (Poly):
                Secret key.
            C (int):
                Number of coefficients to bootstrap. Must be a power of two and
                at least 2.
            g (int):
                Grouping parameter (logarithm of radix) for SlotToCoeff.
            print_messages (bool, optional):
                Whether to print information messages. Defaults to False.
        """
        if 2 * C * cls.h > cls.N:
            raise ValueError("2 * C * h is bounded by N.")
        if C & (C - 1) != 0:
            raise ValueError("C must be a power of two.")
        if C < 2:
            raise ValueError("C must be at least 2.")

        cls._check_key_gen()

        # Generate the bootstrapping keys

        if print_messages:
            print("Generating bootstrapping keys...")
        cs_list = []
        cls.C = C
        B = cls.N // cls.h
        log_C_half = log(C // 2, 2)
        sk_coeffs = cls.sk.coeffs
        for u in range(2 * cls.C):
            sk_tilde = np.zeros(cls.N // 2, dtype=np.complex128)
            for k in range(B // (2 * cls.C)):
                for b in range(cls.h):
                    for a in range(cls.C):
                        sk_tilde[k * cls.h * cls.C + b * cls.C + a] = int(
                            sk_coeffs[b * B + u * B // (2 * cls.C) + k]
                        )
            sk_tilde = ext_bit_rev_vector(sk_tilde, cls.C // 2, log_C_half)
            cs_list.append(
                cls.enc_poly_with_sk(
                    cls.encode(sk_tilde, is_boot=True), sk, is_boot=True
                )
            )
        cls.cs_list = cs_list

        # Generate and encode vector eta

        if print_messages:
            print("Encoding relevant vector as polynomial...")
        eta = np.zeros(cls.C, np.complex128)
        for i in range(cls.C):
            eta[i] = 1 if i < cls.C // 2 else 1j
        cls.encoded_eta = cls.encode(eta, is_boot=True)

        # Generate and encode matrices for SlotToCoeff

        if print_messages:
            print(
                "Generating and encoding matrices required for SlotToCoeff..."
            )
        g = min(g, log(cls.C, 2) - 1)
        grouped_F_list = get_grouped_F(cls.C // 2, g, False).copy()
        cls.grouped_poly_F_list = [
            cls.get_poly_matrix(A, is_boot=True) for A in grouped_F_list
        ]

        # Generate missing switching keys

        if print_messages:
            print("Generating missing switching keys...")
        cls.get_galois_swk(
            -1, sk, cls.moduli_boot[-2], cls.P
        )  # For conjugation
        rotation_indices = [
            2**i for i in range(log(cls.N, 2))
        ]  # For rotation by powers of two
        for poly_matrix in cls.grouped_poly_F_list:
            rotation_indices += cls.get_BSGS_rotation_indices(poly_matrix)
        for i in rotation_indices:
            cls.get_galois_swk(5**i, sk, cls.moduli_boot[-2], cls.P)

        if print_messages:
            print("The SPRU bootstrapping configuration is done!")

    def SPRU(self):
        """
        Apply the SPRU bootstrapping procedure. The ciphertext self is
        refreshed by increasing its level. Only the coefficients indexed by
        multiples of N / C are bootstrapped.

        Returns:
            CKKS:
                A new ciphertext with increased level.
        """
        # Coefficient encodings

        B = self.N // self.h
        outer_factor = (self.q0 / (4 * self.delta * np.pi)) ** (1 / self.h)
        inner_factor = 2 * np.pi * 1j / self.q0
        a_coeffs = self.a.coeffs
        b_coeffs = self.b.coeffs
        E_list = []
        for u in range(2 * self.C):
            e = np.zeros(self.N // 2, dtype=np.complex128)
            for k in range(B // (2 * self.C)):
                for b in range(self.h):
                    for a in range(self.C):
                        i = b * B + u * B // (2 * self.C) + k
                        if self.N // self.C * a >= i:
                            e_coeff = a_coeffs[self.N // self.C * a - i]
                        else:
                            e_coeff = -a_coeffs[
                                self.N + self.N // self.C * a - i
                            ]
                        if i == 0:
                            e_coeff += b_coeffs[self.N // self.C * a]
                        e[k * self.h * self.C + b * self.C + a] = int(e_coeff)
            e = np.exp(e * inner_factor) * outer_factor
            e = ext_bit_rev_vector(e, self.C // 2)
            E_list.append(self.encode(e, is_boot=True))

        # Homomorphic operations until SlotToCoeff

        ct = self.enc_poly_without_error(0, is_boot=True)
        for u in range(2 * self.C):
            ct += self.cs_list[u] @ E_list[u]
        ct = ct.trace(self.N // 2, self.h * self.C)
        ct = ct.product(self.h * self.C, self.C)
        ct = ct - ct.conjugate()
        ct = -(1j * ct)

        # SlotToCoeff

        ct = self.encoded_eta @ ct
        ct = ct.trace(self.C, self.C // 2)
        for A in self.grouped_poly_F_list:
            ct = ct.BSGS_left_mult(A)
        ct = ct.boot_to_nonboot()

        return ct
