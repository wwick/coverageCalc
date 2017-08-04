import os
from subprocess import call
import numpy as np
import pandas as pd
import time
from multiprocessing import Pool

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

filePath = "/users/wwick/adrianBAM"
ref = "AE004437.1"
numThreads = 64
numTPs = 4
numReps = 3

genes = {}
genes["VNG_0780H"] = ["586952", "587128"]
genes["VNG_1130H"] = ["858687", "858911"]
genes["VNG_2446H"] = ["1835869", "1836027"]
genes["VNG_0287H"] = ["230608", "230883"]

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

def coverageComputer(gene):
    plt.figure()
    for tp in xrange(1, 5):
        BAMfiles = []
        for rep in xrange(1, 4):
            BAMfiles.append(BAMfile(rep, tp))
        output = outFile(gene, tp)
        stdout = open(output, "w")
        stderr = open(errFile(gene, tp), "w")
        command = ["samtools", "mpileup", "-d 8000", "-r", genePos(gene),
            BAMfiles[0], BAMfiles[1], BAMfiles[2]]
        call(command, stdout = stdout, stderr = stderr)
        stdout.close()
        stderr.close()
        pileup = np.loadtxt(output, skiprows = 0, usecols = [1,3])
        pileup[:,0] = pileup[:,0] - pileup[0,0] + 1
        pileup[:,1] = pileup[:,1] / pileup[:,1].sum() * 1e6
        plt.plot(pileup[:,0], pileup[:,1], label = "Time Point " + str(tp))
    plt.xlabel('relative position')
    plt.ylabel('reads per million')
    plt.title(gene)
    plt.legend()
    plt.savefig(figFile(gene, tp), dpi = 600)

t0 = time.time()
hydra = Pool(numThreads)
hydra.map(coverageComputer, genes)
t1 = time.time()
print(t1 - t0)
