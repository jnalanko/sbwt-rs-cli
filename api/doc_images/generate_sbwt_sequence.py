import sys

seqs = ["TGTTTG", "TTGCTAT", "ACGTAGTATAT", "TGTAAA"]
k = 4

kmers = []
source_kmer_set = set()

for T in seqs:
    for i in range(len(T)-k+1):
        kmers.append(T[i:i+k])

for T in seqs:
    for i in range(len(T)-k+1):
        x = T[i:i+k]
        if all([(c + x[0 : k-1]) not in kmers for c in "ACGT"]):
            source_kmer_set.add(x)
            for j in range(k-1):
                x = '$' + x[0 : k-1]
                kmers.append(x)
kmers.append("$"*k)

kmers = list(set(kmers)) # Remove duplicates
kmers.sort(key = lambda x : x[::-1]) # colex

kmer_set = set(kmers)

import sys

for i, x in enumerate(kmers):
    outlabels = []
    x_other_set = [d + x[1:] for d in "$ACGT" if (d + x[1:]) in kmer_set]
    for c in "ACGT":
        y = (x[1:] + c)
        if y in kmer_set:
            j = kmers.index(y)
            if x == min(x_other_set): 
                outlabels.append(c)
    #sys.stdout.write("| " + x + " | {" + ",".join(outlabels) + "} |\n")
    if i > 0:
        sys.stdout.write(", ")
    sys.stdout.write("{" + ",".join(outlabels) + "}")
print()
            
