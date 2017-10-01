# filename = "rosalind_ba3c.txt"
# file = open(filename, "r")
kmers = []
dna = 'AAGATTCTCTAC'
k = 4
for x in range(0, (len(dna)-k)+1):
    temp = dna[x:(x+k)]
    kmers.append(temp)

# for line in file:
#     kmers.append(line.rstrip('\r\n'))

kmers.sort()





