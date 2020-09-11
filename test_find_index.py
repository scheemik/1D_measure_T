import numpy as np
import helper_functions as hf
import pdb

test_array = np.linspace(-1, 1, 1024)

z_array = np.genfromtxt('z_array.csv', delimiter=',')

print('z_array.dtype',z_array.dtype)
print(z_array)
print('test_array.dtype',test_array.dtype)

test_this = 1.464839064821421566e-03

print(hf.find_nearest_index(np.flip(z_array), 0.002, 0.1))

# pdb.run('print(hf.find_nearest_index(z_array, test_this, 0.1))')
