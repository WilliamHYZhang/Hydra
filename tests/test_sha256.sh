#!/usr/bin/env bash

cd ../
echo "compiling rust polycommit..."
cd deps/poly-commit/ && cargo build --release
echo "done."
cd ../../
echo "compiling subcircuit..."
cd protocols/sha256/ && g++ -fopenmp -o subcircuit sha256_subcircuit.cpp
echo "done."
cd ../../
echo "compiling hydra..."
cd protocols/sha256/ && g++ -fopenmp -o hydra sha256_hydra.cpp

echo "done."

echo "running subcircuit tests..."
echo "5383 0" | ./subcircuit > logs/subcircuit_result1.log
echo "68 1" | ./subcircuit > logs/subcircuit_result2.log
echo "10 1" | ./subcircuit > logs/subcircuit_result2.log

echo "done."

echo "running hydra tests..."
echo "0 0" | ./hydra > logs/hydra_result1.log
echo "0 1" | ./hydra > logs/hydra_result2.log
echo "done."