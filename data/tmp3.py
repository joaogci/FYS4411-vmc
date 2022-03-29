from process_data import process_data
import os
import sys

directory = "./part_b_sl_comp/"
folders = os.listdir(directory)
try:
    folders.remove(".DS_Store")
except:
    ...
folders.sort()

sl = [0.001, 0.01, 0.1, 1.0, 10.0]

for i, folder in enumerate(folders):
    process_data(directory + folder + "/", "./data_sl" + str(sl[i]) + ".csv")
                
