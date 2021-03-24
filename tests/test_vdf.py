#Verifiable Delay Function Tests
import os
import random 
import hashlib

print('compiling rust polycommit...')
os.chdir('../deps/poly-commit/')
os.system('cargo build --release')
print('done.')

print('compiling subcircuit...')
os.chdir('../../protocols/sha256/')
os.system('g++ -fopenmp -o subcircuit sha256_subcircuit.cpp')
print('done.')

print('compiling hydra...')
os.system('g++ -fopenmp -o hydra sha256_hydra.cpp')
print('done.')

def test(n):
  start = random.getrandbits(256).to_bytes(256, byteorder='big')
  hash = int(hashlib.sha256(start).hexdigest(), 16)
  for i in range(n):
    print(f'testing {n} repeats, iteration {i}')
    payload = bin(hash)[2:]
    payload = '0'*(512 - len(payload)) + payload
 
    print('hydra...')
    os.system(f'echo "0 1 {payload}" | ./hydra > logs/hydra_result_{i}_{n}.log')
    os.system(f'echo "0 0 {payload}" | ./hydra > logs/hydra_control_result_{i}_{n}.log')

    print('subcircuit...')
    os.system(f'echo "68 1 {payload}" | ./subcircuit > logs/subcircuit_result_{i}_{n}.log')
    os.system(f'echo "5383 0 {payload}" | ./subcircuit > logs/subcircuit_control_result_{i}_{n}.log')

    hash = int(hashlib.sha256(hash.to_bytes(256, byteorder='big')).hexdigest(), 16)

test(25)
