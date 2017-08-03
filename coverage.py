import os
from subprocess import call
import numpy as np
import pandas as pd
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

filePath = "/users/wwick/adrianBAM"
ref = "AE004437.1"

def BAMfile(rep, tp):
    return (filePath + "/BAM/rbf.rep." + str(rep) + ".tp." + str(tp) +
        "/Aligned.sortedByCoord.out.bam")

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

genes = {}
genes["VNG_0780H"] = ["586952", "587128"]
genes["VNG_1130H"] = ["858687", "858911"]
genes["VNG_2446H"] = ["1835869", "1836027"]
genes["VNG_0287H"] = ["230608", "230883"]

for gene in genes:
    plt.figure()
    for tp in range(1, 5):
        BAMfiles = []
        for rep in range(1, 4):
            BAMfiles.append(BAMfile(rep, tp))
        output = outFile(gene, tp)
        stdout = open(output, "w")
        stderr = open(errFile(gene, tp), "w")
        command = ["samtools", "mpileup", "-d 8000", "-r", genePos(gene),
            BAMfiles[0], BAMfiles[1], BAMfiles[2]]
        t0 = time.time()
        call(command, stdout = stdout, stderr = stderr)
        t1 = time.time()
        print(t1 - t0)
        stdout.close()
        stderr.close()
        pileup = np.loadtxt(output, skiprows = 0, usecols = [1,3])
        x = pileup[:,0]
        x = x - x[0] + 1
        y = pileup[:,1]
        plt.plot(x,y)
    plt.xlabel('relative position')
    plt.ylabel('read depth')
    plt.title(gene)
    plt.savefig(figFile(gene, tp))
