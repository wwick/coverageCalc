import os
from subprocess import call
import numpy as np
import time
from multiprocessing import Process

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

filePath = "/users/wwick/adrianBAM"
ref = "AE004437.1"
numTP = 4
numReps = 3

genes = {}

# Outliers
genes["VNG_1130H"] = ["858687", "858911"] # all
genes["VNG_2446H"] = ["1835869", "1836027"] # all
genes["VNG_0287H"] = ["230608", "230883"] # 1,2,3
genes["VNG_1207C"] = ["907027", "907659"] # 3,4

# Controls
genes["VNG_1727G"] = ["1277691", "1278269"] # 1
genes["VNG_1617H"] = ["1205478", "1205966"] # 2
genes["VNG_0713C"] = ["537212", "537907"] # 3
genes["VNG_2140G"] = ["1575479", "1575841"] # 4

def BAMfile(rep, tp):
    return (filePath + "/BAM/rbf.rep." + str(rep) + ".tp." + str(tp) +
        "/Aligned.sortedByCoord.out.bam")

def tempFile(rep, tp):
    path = filePath + "/temp/"
    if not os.path.exists(path):
        os.makedirs(path)
    return (path + "rep." + str(rep) + ".tp." + str(tp) + ".temp.txt")

def outFile(gene, tp):
    path = filePath + "/output/" + gene
    if not os.path.exists(path):
        os.makedirs(path)
    return (path + "/tp." + str(tp) + ".out.txt")

def errFile(gene, tp):
    path = filePath + "/error/" + gene
    if not os.path.exists(path):
        os.makedirs(path)
    return (path + "/tp." + str(tp) + ".err.txt")

def figFile(gene, tp):
    path = filePath + "/figures/"
    if not os.path.exists(path):
        os.makedirs(path)
    return (path + gene + ".plot.png")

def genePos(gene):
    return (ref + ":" + genes[gene][0] + "-" + genes[gene][1])

def readCalc(rep, tp):
    input = BAMfile(rep, tp)
    stdout = open(tempFile(rep, tp), "w")
    command = ["samtools", "view", "-c", "-F", "260", input]
    call(command, stdout = stdout)

def coverageCalc(gene):
    plt.figure()
    for tp in range(1, numTP + 1):
        inputs = []
        for rep in range(1, numReps + 1):
            inputs.append(BAMfile(rep, tp))
        output = outFile(gene, tp)
        stdout = open(output, "w")
        stderr = open(errFile(gene, tp), "w")
        command = ["samtools", "depth", "-r", genePos(gene),
            inputs[0], inputs[1], inputs[2]]
        call(command, stdout = stdout, stderr = stderr)
        stdout.close()
        stderr.close()
        coverage = np.loadtxt(output, skiprows = 0, usecols = [1,2])
        x = coverage[:,0]
        y = coverage[:,1]
        x = x - x[0] + 1
        y = y / numReads[tp - 1] * 1e6
        # y = y / np.amax(y) * 1e2
        plt.plot(x, y, label = "Time Point " + str(tp))
    plt.xlabel('relative position')
    plt.ylabel('read depth')
    plt.title(gene)
    plt.legend()
    plt.savefig(figFile(gene, tp), dpi = 600)


t0 = time.time()

tasks = []
# for tp in range(1, numTP + 1):
#      for rep in range(1, numReps + 1):
#             task = Process(target = readCalc, args = (rep, tp))
#             tasks.append(task)
#             task.start()
# for task in tasks:
#     task.join()

numReads = []
for tp in range(1, numTP + 1):
    values = []
    for rep in range(1, numReps + 1):
        path = tempFile(rep, tp)
        file = open(path)
        values.append(int(file.read()))
        file.close()
    numReads.append(np.mean(values))

del tasks [:]
for gene in genes:
    task = Process(target = coverageCalc, args = (gene,))
    tasks.append(task)
    task.start()
for task in tasks:
    task.join()

t1 = time.time()
print(t1 - t0)
