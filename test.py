print "hellp world"

string2 = "ATGTAAATACCGGGTCCGAGGTGCGTCCAGTGGCCGGCGGAGCTGATTCCCTAACTATTCTGTATAAGCCTCCTGCCGCTGTAGTTATCTACGAACTCTACTACAAATGGTCAGGGTCCTATGTAGCGTAAACGACTTCGTCCTCTCTTCCCGTGACCGTTGAGCTTGGTGGAACCCGCACGAGATTTATGCACGTGTGA"
k = 7
A = [0.179, 0.179, 0.286, 0.214, 0.286, 0.286, 0.214]
C = [0.321, 0.286, 0.286, 0.393, 0.214, 0.25, 0.357]
G = [0.25, 0.321, 0.25, 0.179, 0.25, 0.357, 0.143]
T = [0.25, 0.214, 0.179, 0.214, 0.25, 0.107, 0.286]
string2len = len(string2)
probability = 0
outString = ""
position = 0
while position < string2len:
    testString = string2[position:(position + k)]
    pos = 0
    sumOfprob = 0
    if len(testString) < k:
        break
    for char in testString:
        if char == "A":
            sumOfprob += A[pos]
        if char == "G":
            sumOfprob += G[pos]
        if char == "T":
            sumOfprob += T[pos]
        if char == "C":
            sumOfprob += C[pos]
        pos += 1
    pos = 0
    if sumOfprob > probability:
        probability = sumOfprob
        outString = testString
    position += 1

print outString


