import os
from subprocess import call
import numpy as np
import time
from multiprocessing import Process
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

file_path = "/users/wwick/adrianBAM"
temp_path = file_path + "/temp"
results_path = file_path + "/results"
data_path = file_path + "/data"

fasta = data_path + "/halo.fasta"
meme_fasta = temp_path + "/meme.fasta"
meme_out = results_path + "/meme_out/"

ref = "AE004437.1"
num_tp = 4
num_reps = 3

genes = {}

# Outliers
genes["VNG_1130H"] = ["858687", "858911"] # all
genes["VNG_2446H"] = ["1835869", "1836027"] # all
genes["VNG_0287H"] = ["230608", "230883"] # 1,2,3
genes["VNG_1207C"] = ["907027", "907659"] # 3,4
# genes["VNG_6322H"] = ["254500", "254730"] # 4
# genes["VNG_0019H"] = ["15310", "15408"] # 3,4


# Controls
# genes["VNG_1727G"] = ["1277691", "1278269"] # 1
# genes["VNG_2140G"] = ["1575479", "1575841"] # 4

def bam_file(rep, tp):
    return (data_path + "/rbf.rep." + str(rep) + ".tp." + str(tp) +
        "/Aligned.sortedByCoord.out.bam")

def num_readsFile(rep, tp):
    path = temp_path + "/numReads/"
    if not os.path.exists(path):
        os.makedirs(path)
    return (path + "rep." + str(rep) + ".tp." + str(tp) + ".temp.txt")

def pileup_file(gene, tp):
    path = temp_path + "/pileups/" + gene + "/"
    if not os.path.exists(path):
        os.makedirs(path)
    return (path + "tp." + str(tp) + ".out.txt")

def err_file(gene, tp):
    path = temp_path + "/error/" + gene + "/"
    if not os.path.exists(path):
        os.makedirs(path)
    return (path + "tp." + str(tp) + ".err.txt")

def fig_file(gene):
    path = results_path + "/figures/"
    if not os.path.exists(path):
        os.makedirs(path)
    return (path + gene + ".plot.png")

def fig_fileTP(gene, tp):
    path = results_path + "/figures/Reps/TP" + str(tp) + "/"
    if not os.path.exists(path):
        os.makedirs(path)
    return (path + gene + ".plot.png")

def gene_pos(gene):
    return (ref + ":" + genes[gene][0] + "-" + genes[gene][1])

def readCalc(rep, tp):
    input = bam_file(rep, tp)
    stdout = open(num_readsFile(rep, tp), "w")
    command = ["samtools", "view", "-c", input]
    call(command, stdout = stdout)

def coverage_calc(gene, tp):
    bams = []
    for rep in range(1, num_reps + 1):
        bams.append(bam_file(rep, tp))
    output = pileup_file(gene, tp)
    stdout = open(output, "w")
    stderr = open(err_file(gene, tp), "w")
    command = ["samtools", "mpileup", "-d 8000", "-r", gene_pos(gene), "-f",
        fasta, bams[0], bams[1], bams[2]]
    call(command, stdout = stdout, stderr = stderr)
    stdout.close()
    stderr.close()

def graph_tp(gene, tp):
    plt.figure()
    data = pileup_file(gene, tp)
    coverage = np.loadtxt(data, skiprows = 0, usecols = [1,3,6,9])
    x = coverage[:,0]
    x = x - x[0] + 1
    for rep in range(1, num_reps + 1):
        y = coverage[:,rep]
        scaling_factor = num_reads[tp - 1][rep - 1] / 1e6
        y /= scaling_factor
        plt.plot(x, y, label = "Rep " + str(rep))
    plt.xlabel('Relative Position')
    plt.ylabel('Normalized Reads**')
    plt.title(gene + " TP " + str(tp))
    plt.legend()
    plt.savefig(fig_fileTP(gene, tp), dpi = 600)

def graph(gene):
    plt.figure()
    ax = plt.subplot()
    gene_df = pd.DataFrame()
    for tp in range(1, num_tp + 1):
        df = pd.read_csv(pileup_file(gene, tp), sep = '\t', header = None,
            names = ['loci', 'base', 'rep 1', 'rep 2', 'rep 3'],
            usecols = [1, 2, 3, 6, 9])
        for rep in range(1, num_reps + 1):
            name = 'rep ' + str(rep)
            scaling_factor = num_reads[tp-1][rep-1] / 1e6
            df[[name]] /= scaling_factor
        name = 'tp ' + str(tp)
        gene_df[name] = df[['rep 1', 'rep 2', 'rep 3']].mean(axis = 1)
    gene_df['relative position'] = gene_df.index + 1
    for tp in range(1, num_tp + 1):
        gene_df.plot('relative position', 'tp ' + str(tp), ax = ax)
    ax.legend(['TP 1', 'TP 2', 'TP 3', 'TP 4'])
    ax.set_title(gene)
    ax.set_xlabel('Relative Position')
    ax.set_ylabel('Reads Per Million Reads in BAM')
    plt.savefig(fig_file(gene), dpi = 600)

def meme():
    meme_in = open(meme_fasta, 'w')
    for gene in genes:
        gene_df = pd.DataFrame()
        for tp in range(1, num_tp + 1):
            tp_df = pd.read_csv(pileup_file(gene,tp), sep = '\t', header = None,
                names = ['loci', 'base', 'rep 1', 'rep 2', 'rep 3'],
                usecols = [1, 2, 3, 6, 9])
            for rep in range(1, num_reps + 1):
                name = 'rep ' + str(rep)
                scaling_factor = num_reads[tp-1][rep-1]
                tp_df[[name]] /= scaling_factor
            # tp_df['mean'] = tp_df[['rep 1', 'rep 2', 'rep 3']].sum(axis = 1)
            # criteria = tp_df['mean'].max() / 2
            # index = tp_df.ix[tp_df['mean'] > criteria].index
            # reach = 0
            # start = min(index) - reach
            # end = max(index) + reach
            name = 'tp ' + str(tp)
            # gene_df[name] = tp_df.loc[range(start,end),'base']
            gene_df[name] = tp_df['base']
        title_line = '>' + gene + '\n'
        sequence_line = ''.join(gene_df['tp 1'].values) + '\n'
        meme_in.write(title_line)
        meme_in.write(sequence_line)
    meme_in.close()
    # command = ['meme', meme_fasta, '-dna', '-evt', str(5e-2), '-oc', meme_out]
    command = ['meme', meme_fasta, '-dna', '-oc', meme_out]
    call(command)

t0 = time.time()
tasks = []

# for tp in range(1, num_tp + 1):
#      for rep in range(1, num_reps + 1):
#             task = Process(target = readCalc, args = (rep, tp))
#             tasks.append(task)
#             task.start()
# for task in tasks:
#     task.join()
# del tasks [:]

num_reads = []
for tp in range(1, num_tp + 1):
    values = []
    for rep in range(1, num_reps + 1):
        path = num_readsFile(rep, tp)
        file = open(path)
        values.append(int(file.read()))
        file.close()
    num_reads.append(values)

for gene in genes:
    for tp in range(1, num_tp + 1):
        task = Process(target = coverage_calc, args = (gene, tp))
        tasks.append(task)
        task.start()
for task in tasks:
    task.join()
del tasks [:]

for gene in genes:
    task = Process(target = graph, args = (gene,))
    tasks.append(task)
    task.start()
    # for tp in range(1, num_tp + 1):
    #     task = Process(target = graph_tp, args = (gene, tp))
    #     tasks.append(task)
    #     task.start()
for task in tasks:
    task.join()

meme()

t1 = time.time()
print(t1 - t0)
