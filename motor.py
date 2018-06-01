import numpy as np
import string

def get_data(filename):
    z = open('motor_curves/'+filename+'.txt')
    time_force=[]
    for line in z:
        try:
            float(line[0])
            time_force.append(map(float, line.split()))
        except ValueError:
            continue
    z.close()
    return time_force