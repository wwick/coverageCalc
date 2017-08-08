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
fasta = filePath + "/data/halo.fasta"
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
    return (filePath + "/data/rbf.rep." + str(rep) + ".tp." + str(tp) +
        "/Aligned.sortedByCoord.out.bam")

def tempFile(rep, tp):
    path = filePath + "/temp/"
    if not os.path.exists(path):
        os.makedirs(path)
    return (path + "rep." + str(rep) + ".tp." + str(tp) + ".temp.txt")

def outFile(gene, tp):
    path = filePath + "/output/" + gene + "/"
    if not os.path.exists(path):
        os.makedirs(path)
    return (path + "tp." + str(tp) + ".out.txt")

def errFile(gene, tp):
    path = filePath + "/error/" + gene + "/"
    if not os.path.exists(path):
        os.makedirs(path)
    return (path + "tp." + str(tp) + ".err.txt")

def figFile(gene):
    path = filePath + "/figures/"
    if not os.path.exists(path):
        os.makedirs(path)
    return (path + gene + ".plot.png")

def figFileTP(gene, tp):
    path = filePath + "/figures/TP" + str(tp) + "/"
    if not os.path.exists(path):
        os.makedirs(path)
    return (path + gene + ".plot.png")

def genePos(gene):
    return (ref + ":" + genes[gene][0] + "-" + genes[gene][1])

def readCalc(rep, tp):
    input = BAMfile(rep, tp)
    stdout = open(tempFile(rep, tp), "w")
    command = ["samtools", "view", "-c", input]
    call(command, stdout = stdout)

def coverageCalc(gene, tp):
    bams = []
    for rep in range(1, numReps + 1):
        bams.append(BAMfile(rep, tp))
    output = outFile(gene, tp)
    stdout = open(output, "w")
    stderr = open(errFile(gene, tp), "w")
    command = ["samtools", "mpileup", "-d 8000", "-r", genePos(gene), "-f",
        fasta, bams[0], bams[1], bams[2]]
    call(command, stdout = stdout, stderr = stderr)
    stdout.close()
    stderr.close()

def graphTP(gene, tp):
    plt.figure()
    data = outFile(gene, tp)
    coverage = np.loadtxt(data, skiprows = 0, usecols = [1,3,6,9])
    x = coverage[:,0]
    x = x - x[0] + 1
    lineStyles = ['dashed', 'dashdot', 'dotted']
    for rep in range(1, numReps + 1):
        y = coverage[:,rep]
        y = y / numReads[tp - 1][rep - 1] * 1e6
        plt.plot(x, y, label = "Rep " + str(rep))
    plt.xlabel('relative position')
    plt.ylabel('reads per mapped reads in .bam times 1e6')
    plt.title(gene + " Time Point " + str(tp))
    plt.legend()
    plt.savefig(figFileTP(gene, tp), dpi = 600)


# def coverageCalc(gene):
#     plt.figure()
#     for tp in range(1, numTP + 1):
#         inputs = []
#         for rep in range(1, numReps + 1):
#             inputs.append(BAMfile(rep, tp))
#         output = outFile(gene, tp)
#         stdout = open(output, "w")
#         stderr = open(errFile(gene, tp), "w")
#         command = ["samtools", "mpileup", "-r", genePos(gene), "-f", fasta,
#             # inputs[0], inputs[1], inputs[2]]
#             inputs[0]]
#         call(command, stdout = stdout, stderr = stderr)
#         stdout.close()
#         stderr.close()
#         coverage = np.loadtxt(output, skiprows = 0, usecols = [1,3])
#         x = coverage[:,0]
#         y = coverage[:,1]
#         x = x - x[0] + 1
#         y = y / numReads[tp - 1] * 1e6
#         # y = y / y.sum()
#         # print(len(x), tp , gene)
#         # y = y / np.amax(y) * 1e2
#         plt.plot(x, y, label = "Time Point " + str(tp))
#     plt.xlabel('relative position')
#     plt.ylabel('read depth')
#     plt.title(gene)
#     plt.legend()
#     plt.savefig(figFile(gene), dpi = 600)


t0 = time.time()
tasks = []

for tp in range(1, numTP + 1):
     for rep in range(1, numReps + 1):
            task = Process(target = readCalc, args = (rep, tp))
            tasks.append(task)
            task.start()
for task in tasks:
    task.join()
del tasks [:]

numReads = []
for tp in range(1, numTP + 1):
    values = []
    for rep in range(1, numReps + 1):
        path = tempFile(rep, tp)
        file = open(path)
        values.append(int(file.read()))
        file.close()
    numReads.append(values)

for gene in genes:
    for tp in range(1, numTP + 1):
        task = Process(target = coverageCalc, args = (gene, tp))
        tasks.append(task)
        task.start()
for task in tasks:
    task.join()
del tasks [:]

for gene in genes:
    for tp in range(1, numTP + 1):
        task = Process(target = graphTP, args = (gene, tp))
        tasks.append(task)
        task.start()
for task in tasks:
    task.join()

t1 = time.time()
print(t1 - t0)
