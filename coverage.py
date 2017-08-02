import os
from subprocess import call

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

def genePos(gene):
    return (ref + ":" + genes[gene][0] + "-" + genes[gene][1])

genes = {}
genes["VNG_0780H"] = ["586952", "587128"]
genes["VNG_1130H"] = ["858687", "858911"]
genes["VNG_2446H"] = ["1835869", "1836027"]
genes["VNG_0287H"] = ["230608", "230883"]

for gene in genes:
    for tp in range(1,5):
        BAMfiles = []
        for rep in range(1,4):
            BAMfiles.append(BAMfile(rep, tp))
        stdout = open(outFile(gene, tp), "w")
        stderr = open(errFile(gene, tp), "w")
        command = ["samtools", "mpileup", "-d 8000", "-r", genePos(gene),
            BAMfiles[0], BAMfiles[1], BAMfiles[2]]
        call(command, stdout = stdout, stderr = stderr)
        stdout.close()
        stderr.close()
