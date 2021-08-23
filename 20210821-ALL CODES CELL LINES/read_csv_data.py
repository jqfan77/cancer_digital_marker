# Last updated on Aug 21th, 2021,  fjq19@mails.tsinghua.edu.cn



import numpy as np
import csv


# load csv data
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

# load csv data, auto detect the delimiter
def read_csv_data_autodelimiter(path, row=-1,firstline=True):
    
    delimiter_list = [',', ' ', '\t']
    
    f = open(path, 'r')
    
    for the_delimiter in delimiter_list:
        
        flag_ok = False
        
        t=[]
        
        f = open(path, 'r')
        d = csv.reader(f, delimiter=the_delimiter)
        
        for line in d:
            if len(line)>=2:
                flag_ok = True
            t.append(line)
        f.close()
        
        if not flag_ok:
            continue
        
        if firstline==False:
            del t[0]
        if row==-1:
            return t
        else:
            return t[row]
    
    return read_csv_data(path, row, firstline)
