

filename = "rosalind_ba3c.txt"
file = open(filename, "r")
kmers = []

for line in file:
    kmers.append(line.rstrip('\r\n'))

kmers.sort()
k = len(kmers[0])
dna = []
prefix_and_sufix = {}
kmer_edges = {}
for kmer in kmers:
    preffix = kmer[0:k-1]
    suffix = kmer[len(kmer)-k+1: len(kmer)]
    up_date = {kmer: [preffix, suffix]}
    prefix_and_sufix.update(up_date)
for entry in prefix_and_sufix:
    pre = prefix_and_sufix[entry][0]
    suff = prefix_and_sufix[entry][1]
    edges = []
    for e in prefix_and_sufix:
        if suff == prefix_and_sufix[e][0]:
            edges.append(e)
    temp = {entry: edges}
    kmer_edges.update(temp)

for x in kmers:
    tempm = x
    tempm += " -> "
    for z in kmer_edges[x]:
        tempm += z
        tempm += ','
    print(tempm)

