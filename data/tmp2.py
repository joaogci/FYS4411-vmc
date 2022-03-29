from process_data import process_data
import os
import sys

directory = "./part_b_test/"
folders = os.listdir(directory)
try:
    folders.remove(".DS_Store")
except:
    ...
folders.sort()

N = [100, 100, 100, 10, 10, 10, 1, 1, 1, 500, 500, 500]
d = [1, 2, 3]

for i, folder in enumerate(folders):
    process_data(directory + folder + "/", "./data_N" + str(N[i]) + "_d" + str(d[i % 3]) + ".csv")




