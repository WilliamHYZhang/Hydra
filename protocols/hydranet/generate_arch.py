with open('arch100.txt', 'w') as f1:
    f1.write('512\n') #batch size
    f1.write('256\n') #input size
    for i in range(100):
        f1.write('512\n') #hidden layers
    f1.write('128\n') #output layer

with open('arch1000.txt', 'w') as f1:
    f1.write('512\n') #batch size
    f1.write('256\n') #input size
    for i in range(1000):
        f1.write('512\n') #hidden layers
    f1.write('128\n') #output layer