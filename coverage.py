import numpy as np
from subprocess import call

filePath = "/users/wwick/adrianBAM/"
BAMfile = filePath + "BAM/rbf.rep.1.tp.1/Aligned.sortedByCoord.out.bam"
outFile = filePath + "output/out.txt"

stdout = open(outFile, "w")

pos = "AE004437.1:586952-587128"
command = ["samtools", "mpileup", "-r", pos, BAMfile]
print(command)
call(command, stdout = stdout)

stdout.close()
