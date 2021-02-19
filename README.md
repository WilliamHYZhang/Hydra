# Hydra

The Hydra Project

## Requirements

- C++11 (g++)
- Rust (cargo)

## Usage

- In `deps/poly-commit/`: compile with `cargo build --release`
- In `protocols/main/protocol_name.cpp`: compile with `g++ -fopenmp -o protocol_name protocol_name.cpp`, run with `./protocol_name`

To run tests (compile, run, and log automatically): `sh test.sh` in `tests/`
