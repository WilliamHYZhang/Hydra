#!/usr/bin/env bash

cd ../
echo "compiling rust polycommit..."
cd deps/poly-commit/ && cargo build --release
echo "done."
cd ../../
echo "compiling subcircuit..."
cd protocols/main/ && g++ -fopenmp -o subcircuit subcircuit.cpp
echo "done."
cd ../../
echo "compiling hydra..."
cd protocols/main/ && g++ -fopenmp -o hydra hydra.cpp

echo "done."

echo "sanity check... (should not hang)"
echo "10 4 10 0" | ./subcircuit
echo "10 4 5 1" | ./subcircuit
echo "10 4 10 0" | ./hydra
echo "10 4 10 1" | ./hydra
echo "done."

echo "testing the final hydra tests..."

echo "1048576 256 0 0" | ./hydra > logs/hydra_result35.log
echo "1048576 256 0 1" | ./hydra > logs/hydra_result36.log

echo "done."