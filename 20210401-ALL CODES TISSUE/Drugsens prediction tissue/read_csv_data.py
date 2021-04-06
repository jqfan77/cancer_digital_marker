import numpy as np
import csv

def read_csv_data(path, row=-1,firstline=True):
    f = open(path, 'r')
    t=[]
    d = csv.reader(f)
    for line in d:
        t.append(line)
    f.close()
    if firstline==False:
        del t[0]
    if row==-1:
        return t
    else:
        return t[row]