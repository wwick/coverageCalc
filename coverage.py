import os
from subprocess import call
import numpy as np
import time
from multiprocessing import Process, Queue
import pandas as pd
import HTSeq
from collections import Counter

import rpy2.robjects as robjects
import rpy2.robjects.packages as rpacks

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
gtf = data_path + "/Halobacterium_salinarum.gtf"
gff = data_path +"/annotation.gff3"
htseq_out = temp_path + "/htseq_out.csv"
deseq_out = temp_path + "/deseq_out.csv"

ref = "AE004437.1"
num_tp = 4
num_reps = 3

genes = {}

# Outliers
genes["VNG_1130H"] = ["858687", "858911"] # all
genes["VNG_2446H"] = ["1835869", "1836027"] # all
genes["VNG_0287H"] = ["230608", "230883"] # 1,2,3
genes["VNG_1207C"] = ["907027", "907659"] # 3,4

# Controls
genes["VNG_1727G"] = ["1277691", "1278269"] # 1
genes["VNG_2140G"] = ["1575479", "1575841"] # 4

def bam_file(rep, tp):
    return (data_path + "/rbf.rep." + str(rep) + ".tp." + str(tp) +
        "/Aligned.sortedByCoord.out.bam")

def num_reads_file(rep, tp):
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

def fig_file_reps(gene):
    path = results_path + "/figures/Reps/"
    if not os.path.exists(path):
        os.makedirs(path)
    return (path + gene + ".plot.png")

def gene_pos(gene):
    return (ref + ":" + genes[gene][0] + "-" + genes[gene][1])

def read_calc(rep, tp):
    input = bam_file(rep, tp)
    stdout = open(num_reads_file(rep, tp), "w")
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

def graph_reps(gene):
    plt.figure()
    ax = plt.subplot()
    colors = ['blue', 'orange', 'green', 'red']
    alpha = 0.4
    linewidth = 2.2

    gene_df = pd.DataFrame()
    for tp in range(1, num_tp + 1):
        df = pd.read_csv(pileup_file(gene, tp), sep = '\t', header = None,
            names = ['rep 1', 'rep 2', 'rep 3'],
            usecols = [3, 6, 9])
        df['relative position'] = df.index + 1
        for rep in range(1, num_reps + 1):
            name = 'rep ' + str(rep)
            scaling_factor = num_reads[(rep-1) + (tp-1) * 3]
            df[[name]] /= scaling_factor
            df.plot('relative position', name, ax = ax, color = colors[tp-1],
                alpha = alpha, label = '_nolegend_')
        name = 'tp ' + str(tp)
        gene_df['TP ' + str(tp)] = df[['rep 1', 'rep 2', 'rep 3']].mean(axis=1)
    gene_df['relative position'] = gene_df.index + 1
    for tp in range(1, num_tp + 1):
        name = 'TP ' + str(tp)
        gene_df.plot('relative position', name, color = colors[tp-1],
            linewidth = linewidth, ax = ax, label = name)
    ax.legend()
    ax.set_title(gene)
    ax.set_xlabel('Relative Position')
    ax.set_ylabel('Reads per Million Total Reads')
    plt.savefig(fig_file_reps(gene), dpi = 600)

def graph(gene):
    plt.figure()
    ax = plt.subplot()
    gene_df = pd.DataFrame()
    for tp in range(1, num_tp + 1):
        df = pd.read_csv(pileup_file(gene, tp), sep = '\t', header = None,
            names = ['rep 1', 'rep 2', 'rep 3'],
            usecols = [3, 6, 9])
        for rep in range(1, num_reps + 1):
            name = 'rep ' + str(rep)
            scaling_factor = num_reads[(rep-1) + (tp-1) * 3]
            df[[name]] /= scaling_factor
        gene_df['TP ' + str(tp)] = df[['rep 1', 'rep 2', 'rep 3']].mean(axis=1)
    gene_df['relative position'] = gene_df.index + 1
    for tp in range(1, num_tp + 1):
        name = 'TP ' + str(tp)
        gene_df.plot('relative position', name, ax = ax, label = name)
    ax.legend()
    ax.set_title(gene)
    ax.set_xlabel('Relative Position')
    ax.set_ylabel('Reads per Million Total Reads')
    plt.savefig(fig_file(gene), dpi = 600)

def meme():
    meme_in = open(meme_fasta, 'w')
    for gene in genes:
        gene_df = pd.DataFrame()
        for tp in range(1, num_tp + 1):
            tp_df = pd.read_csv(pileup_file(gene,tp), sep = '\t', header = None,
                names = ['loci', 'base', 'rep 1', 'rep 2', 'rep 3'],
                usecols = [1, 2, 3, 6, 9])
            tp_df['mean'] = tp_df[['rep 1', 'rep 2', 'rep 3']].mean(axis = 1)
            criteria = tp_df['mean'].max() / 2
            index = tp_df.ix[tp_df['mean'] > criteria].index
            reach = 10
            start = min(index) - reach
            end = max(index) + reach
            name = 'tp ' + str(tp)
            gene_df[name] = tp_df.loc[range(start,end),'base']
        title_line = '>' + gene + '\n'
        print(gene_df)
        sequence_line = ''.join(gene_df['tp 1'].values) + '\n'
        meme_in.write(title_line)
        meme_in.write(sequence_line)
    meme_in.close()
    command = ['meme', meme_fasta, '-dna', '-oc', meme_out]
    call(command)

def htseq():
    gff_file = HTSeq.GFF_Reader(gff, end_included = True)
    genes = HTSeq.GenomicArrayOfSets("auto", stranded = True)

    for gene in gff_file:
        if gene.type == "gene":
            genes[gene.iv] += gene.attr["ID"]

    tasks = []
    results = []
    out_q = Queue()
    for tp in range(1, num_tp + 1):
        for rep in range(1, num_reps + 1):
            task = Process(target = htseq_count, args = (genes, rep, tp, out_q))
            tasks.append(task)
            task.start()
    for i in range(len(tasks)):
        results.append(out_q.get())
    for task in tasks:
        task.join()
    counts = range(len(results))
    for result in results:
        counts[result[1]] = result[0]
    keys_a = set(counts[0].keys())
    for sample in counts:
        keys_b = set(sample.keys())
        keys_a = keys_a & keys_b
    keys = list(keys_a)
    matrix = np.zeros(shape = (len(keys), len(counts)), dtype = np.int_)
    for key_index in range(len(keys)):
        for index in range(len(counts)):
            matrix[key_index, index] = counts[index][keys[key_index]]
    np.savetxt(htseq_out, matrix, delimiter = ",")

def htseq_count(genes, rep, tp, out_q):
    index = (rep-1) + (tp-1) * 3
    count = Counter()
    almnt_file = HTSeq.BAM_Reader(bam_file(rep, tp))
    for almnt in almnt_file:
       gene_ids = set()
       for cigop in almnt.cigar:
           if cigop.type != "M":
              continue
           for iv, val in genes[cigop.ref_iv].steps():
              gene_ids |= val
       if len(gene_ids) == 1:
          gene_id = list(gene_ids)[0]
          count[gene_id] += 1
    output = [count, index]
    out_q.put(output)

def deseq2():

    rpacks.quiet_require('DESeq2')
    robjects.r('countData <- read.csv(\'' + htseq_out + '\')')
    robjects.r('factors = unname(estimateSizeFactorsForMatrix(countData))')
    robjects.r('write.csv(factors, \'' + deseq_out + '\')')

t0 = time.time()
tasks = []

for tp in range(1, num_tp + 1):
     for rep in range(1, num_reps + 1):
            task = Process(target = read_calc, args = (rep, tp))
            tasks.append(task)
            task.start()
for task in tasks:
    task.join()
del tasks [:]

num_reads = []
for tp in range(1, num_tp + 1):
    for rep in range(1, num_reps + 1):
        path = num_reads_file(rep, tp)
        file = open(path)
        num_reads.append(int(file.read())/1e6)
        file.close()

for gene in genes:
    for tp in range(1, num_tp + 1):
        task = Process(target = coverage_calc, args = (gene, tp))
        tasks.append(task)
        task.start()
for task in tasks:
    task.join()
del tasks [:]

htseq()
deseq2()
meme()

factor_df = pd.read_csv(deseq_out)
deseq_factors = list(factor_df['x'])

for gene in genes:
    task = Process(target = graph, args = (gene,))
    task.start()
    tasks.append(task)

    task_reps = Process(target = graph_reps, args = (gene,))
    task_reps.start()
    tasks.append(task_reps)

for task in tasks:
    task.join()

t1 = time.time()
print(t1 - t0)
