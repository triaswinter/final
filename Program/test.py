import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import time
from scipy.optimize import curve_fit
from itertools import count, product
from multiprocessing.dummy import Pool as ThreadPool
import Simulation as SIM
import Matrix as MTX
import mysql.connector

mydb = mysql.connector.connect(
    host="localhost",
    user="root",
    passwd="zoomg21u",
    database="testdb"
)

mycursor = mydb.cursor()

sqlFormula = "INSERT INTO students (name, age) VALUES (%s, %s)"
student1 = ("Tom", 25)

mycursor.execute(sqlFormula, student1)

mydb.commit()
'''def Simulation(para1, para2):
    x=0
    for i in range(1000000):
        x += i*para1+para2
    return x


def run(J, Nt):
    matrix = MTX.create_a(lat_len, para1=J, para2=Nt)
    hmc = SIM.HMC('A', matrix, iteration_len, leapfrog_len, lat_len, Nt)
    configuration = hmc.create_sample_space(0.1)
    return configuration[1592]


my_array1 = (0.1,2,5)
my_array2 = (3,15)
inputpara = product(my_array1, my_array2)
pool = ThreadPool(4)

results = pool.starmap(Simulation, list(inputpara))


lat_len = 3
leapfrog_len = 2
iteration_len = 20000

results = pool.starmap(run, list(inputpara))

print(results)
'''