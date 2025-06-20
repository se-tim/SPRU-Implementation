# SPRU in SageMath

This repository provides a proof-of-concept implementation in SageMath
of the SPRU bootstrapping procedure [1].

## Table of Contents

- [Installation](#installation)
- [Repository Structure](#repository-structure)
- [Usage](#usage)
- [References](#references)

## Installation

1. **Prerequisites**:
Ensure you have SageMath 10.0 or higher installed.
You can download it from the
[SageMath website](https://www.sagemath.org/download.html).
   
2. **Clone the repository**:
   ```bash
   git clone https://github.com/se-tim/SPRU-Implementation.git
   cd SPRU-Implementation
   ```

3. **Initialize the submodules**:
Ensure all submodules are correctly initialized,
as explained here:
[Git submodule inside of a submodule (nested submodules)](https://stackoverflow.com/q/1535524).
Therefore, run the following command inside the cloned repository:
   ```bash
   git submodule update --init --recursive
   ```

## Repository Structure

- **`spru_package/`**:
The main package containing all implementations
related to CKKS and the SPRU procedure:
   - **`ckks_in_sagemath/`**:
   A submodule.
   It implements the main CKKS functionalities,
   including the original bootstrapping procedure.
   Make sure that this is initialized correctly,
   as described in the section [Installation](#installation).
  
   - **`ckks_x.py`**:
   Implements the functionalities required for the SPRU procedure,
   building on the core CKKS operations.

   - **`ext_bit_rev.py`**:
   Implements functions related to extended bit-reversing
   (bit-reversing with only the least significant bits).

- **`test.py`**
 A script demonstrating the SPRU bootstrapping.

## Usage

The `test.py` script demonstrates the SPRU bootstrapping.

1. **Parameter selection**:
Choose default parameters or manually input your own.

2. **Key generation and bootstrapping configuration**:
Generate secret, public and evaluation keys.
Also perform the precomputations required for the SPRU bootstrapping.

3. **Security estimation**:
Estimate the security level of the current parameter set.
By default, the primal hybrid security level is not estimated (to save time).

4. **Encryption**:
Encrypt a random complex vector.

5. **SPRU bootstrapping**:
Refresh a ciphertext at lowest level through bootstrapping.

7. **Precision measurement**:
Evaluate the precision loss introduced during SPRU bootstrapping.

To run the script,
make sure you are in the root directory of the repository,
then execute:
```bash
sage test.py
```
Results will be printed directly to the terminal.

## References

1. Jean-Sébastien Coron and Robin Köstler.
Low-Latency Bootstrapping for CKKS using Roots of Unity, 2025.
https://eprint.iacr.org/2025/651.