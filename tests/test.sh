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

echo "running subcircuit tests..."

#2^16 depth 2^7 width
echo "65536 128 65536 0" | ./subcircuit > logs/subcircuit_result1.log
echo "65536 128 32768 1" | ./subcircuit > logs/subcircuit_result2.log
echo "65536 128 16384 1" | ./subcircuit > logs/subcircuit_result3.log
echo "65536 128 8192 1" | ./subcircuit > logs/subcircuit_result4.log
echo "65536 128 4096 1" | ./subcircuit > logs/subcircuit_result5.log
echo "65536 128 2048 1" | ./subcircuit > logs/subcircuit_result6.log
echo "65536 128 1024 1" | ./subcircuit > logs/subcircuit_result7.log
echo "65536 128 512 1" | ./subcircuit > logs/subcircuit_result8.log
echo "65536 128 256 1" | ./subcircuit > logs/subcircuit_result9.log

#2^16 depth 2^8 width
echo "65536 256 65536 0" | ./subcircuit > logs/subcircuit_result11.log
echo "65536 256 32768 1" | ./subcircuit > logs/subcircuit_result12.log
echo "65536 256 16384 1" | ./subcircuit > logs/subcircuit_result13.log
echo "65536 256 8192 1" | ./subcircuit > logs/subcircuit_result14.log
echo "65536 256 4096 1" | ./subcircuit > logs/subcircuit_result15.log
echo "65536 256 2048 1" | ./subcircuit > logs/subcircuit_result16.log
echo "65536 256 1024 1" | ./subcircuit > logs/subcircuit_result17.log
echo "65536 256 512 1" | ./subcircuit > logs/subcircuit_result18.log
echo "65536 256 256 1" | ./subcircuit > logs/subcircuit_result19.log

#2^20 depth 2^7 width
echo "1048576 128 16777216 0" | ./subcircuit > logs/subcircuit_result21.log
echo "1048576 128 524288 1" | ./subcircuit > logs/subcircuit_result22.log
echo "1048576 128 262144 1" | ./subcircuit > logs/subcircuit_result23.log
echo "1048576 128 131072 1" | ./subcircuit > logs/subcircuit_result24.log
echo "1048576 128 65536 1" | ./subcircuit > logs/subcircuit_result25.log
echo "1048576 128 32768 1" | ./subcircuit > logs/subcircuit_result26.log
echo "1048576 128 16384 1" | ./subcircuit > logs/subcircuit_result27.log
echo "1048576 128 8192 1" | ./subcircuit > logs/subcircuit_result28.log
echo "1048576 128 4096 1" | ./subcircuit > logs/subcircuit_result29.log

#2^20 depth 2^8 width
echo "1048576 256 16777216 0" | ./subcircuit > logs/subcircuit_result31.log
echo "1048576 256 524288 1" | ./subcircuit > logs/subcircuit_result32.log
echo "1048576 256 262144 1" | ./subcircuit > logs/subcircuit_result33.log
echo "1048576 256 131072 1" | ./subcircuit > logs/subcircuit_result34.log
echo "1048576 256 65536 1" | ./subcircuit > logs/subcircuit_result35.log
echo "1048576 256 32768 1" | ./subcircuit > logs/subcircuit_result36.log
echo "1048576 256 16384 1" | ./subcircuit > logs/subcircuit_result37.log
echo "1048576 256 8192 1" | ./subcircuit > logs/subcircuit_result38.log
echo "1048576 256 4096 1" | ./subcircuit > logs/subcircuit_result39.log

echo "done."

echo "running hydra tests..."

#2^16 depth 2^7 width
echo "65536 128 10 0" | ./hydra > logs/hydra_result1.log
echo "65536 128 10 1" | ./hydra > logs/hydra_result2.log
echo "65536 128 20 0" | ./hydra > logs/hydra_result3.log
echo "65536 128 20 1" | ./hydra > logs/hydra_result4.log

#2^16 depth 2^8 width
echo "65536 256 10 0" | ./hydra > logs/hydra_result11.log
echo "65536 256 10 1" | ./hydra > logs/hydra_result12.log
echo "65536 256 20 0" | ./hydra > logs/hydra_result13.log
echo "65536 256 20 1" | ./hydra > logs/hydra_result14.log

#2^20 depth 2^7 width
echo "1048576 128 10 0" | ./hydra > logs/hydra_result21.log
echo "1048576 128 10 1" | ./hydra > logs/hydra_result22.log
echo "1048576 128 20 0" | ./hydra > logs/hydra_result23.log
echo "1048576 128 20 1" | ./hydra > logs/hydra_result24.log

#2^20 depth 2^8 width
echo "1048576 256 10 0" | ./hydra > logs/hydra_result31.log
echo "1048576 256 10 1" | ./hydra > logs/hydra_result32.log
echo "1048576 256 20 0" | ./hydra > logs/hydra_result33.log
echo "1048576 256 20 1" | ./hydra > logs/hydra_result34.log

echo "done."
