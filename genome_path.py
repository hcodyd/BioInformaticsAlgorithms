
filename = "rosalind_ba3b (1).txt"
file = open(filename, "r")
kmers = []

for line in file:
    kmers.append(line.rstrip('\r\n'))


k = 0

dna = ''

dna = kmers[0]

for i in range(1, len(kmers)):
    dna += (kmers[i][len(kmers[i])-1])

print(dna)


