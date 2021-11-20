# Hydra

The Hydra Project (https://www.hydrasec.ml/)

ePrint (extended version): https://eprint.iacr.org/2021/641

Accepted to IEEE TPS 2021.

## Requirements

Linux Environment

- C++11 (g++ version 7.5+)
- Rust (cargo version 1.5+)
- Python (python3 version 3.6.9+)

## Usage

- In `deps/poly-commit/`: compile with `cargo build --release`.
- In `protocols/main/protocol_name.cpp`: compile with `g++ -fopenmp -o protocol_name protocol_name.cpp`, run with `./protocol_name`.

## Testing

To run tests (compile, run, and log automatically): `sh name_of_test.sh` in `tests/`.

VDF tests are run seperately with `python3 test_vdf.py` in `tests/`.
