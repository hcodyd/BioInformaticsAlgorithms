from random import randint
import os



def get_blank_profile(kmer_length, num):
    probability_profile = {'A': [], 'C': [], 'G': [], 'T': []}
    n = 0
    while n < kmer_length:
        probability_profile['A'].insert(n, num)
        probability_profile['C'].insert(n, num)
        probability_profile['G'].insert(n, num)
        probability_profile['T'].insert(n, num)
        n += 1
    return probability_profile


def update_count_profile(k_mer, profile):
    pos = 0
    for char in k_mer:
        profile[char][pos] += 1
        pos += 1


def calculate_prob_profile(profile, k):
    i = 0
    probability_profile = get_blank_profile(k, 1)
    while i < k:
        num_a = profile['A'][i]
        num_g = profile['G'][i]
        num_c = profile['C'][i]
        num_t = profile['T'][i]
        total = num_a + num_c + num_g + num_t
        if num_t != 0:
            probability_profile['T'][i] = num_t/total
        if num_c != 0:
            probability_profile['C'][i] = num_c/total
        if num_g != 0:
            probability_profile['G'][i] = num_g/total
        if num_a != 0:
            probability_profile['A'][i] = num_a/total
        i += 1
    return probability_profile


def test_kmer(kmer, prob_prob):
    probability = 0
    pos = 0
    for char in kmer:
        if pos == 0:
            probability = prob_prob[char][pos]
        else:
            probability = prob_prob[char][pos] * probability
        pos += 1
    return probability


def find_best_motif(dna, count_profile, probability_profile, k):
    first_kmer = dna[0:k]
    best_prob = 0
    best_kmer = ''
    for kmer in range(0, len(dna)):
        kmer_temp = dna[kmer:(kmer + k)]
        if len(kmer_temp) < k:
            break
        kmer_prob = test_kmer(kmer_temp, probability_profile)
        if kmer_prob > best_prob:
            best_prob = kmer_prob
            best_kmer = kmer_temp
    if best_prob == 0:
        best_kmer = first_kmer
    return best_kmer


def get_consensus_sequence(motifs):
    count = get_blank_profile(len(motifs[0]), 0)
    consensus_string = ''
    for i in range(0, len(motifs)):
        update_count_profile(motifs[i], count)
    for j in range(0, len(motifs[0])):
        highest_count = 0
        highest_letter = ''
        if count['A'][j] > highest_count:
            highest_count = count['A'][j]
            highest_letter = 'A'
        if count['C'][j] > highest_count:
            highest_count = count['C'][j]
            highest_letter = 'C'
        if count['T'][j] > highest_count:
            highest_count = count['T'][j]
            highest_letter = 'T'
        if count['G'][j] > highest_count:
            highest_count = count['G'][j]
            highest_letter = 'G'
        consensus_string += highest_letter
    return consensus_string


def hamming_score(motifs):
    motifs_consensus_sequence = get_consensus_sequence(motifs)
    score = 0
    for l in range(0, len(motifs)):
        temp = motifs[l]
        for pos in range(0, len(temp)):
            if temp[pos] != motifs_consensus_sequence[pos]:
                score += 1
    return score

def random_get(dna,k):
    i = 0
    while i < 1:
        rand_spot = randint(0, len(dna))
        temp = dna[rand_spot: (rand_spot + k)]
        if len(temp) == k:
            return temp

def randomized_motif_search(dna, k, t):
    motifs = []
    ii = 0
    for ii in range(0, t):
        motifs.append(random_get(dna[ii], k))
    bestMotifs = motifs

    x = 0
    while x < 1001:
        motifs = []
        count_profile = get_blank_profile(k, 1)
        for motif in bestMotifs:
            update_count_profile(motif, count_profile)
        probability_profile = calculate_prob_profile(count_profile, k)
        for tt in range(0, t):
            temp = find_best_motif(dna[tt], count_profile, probability_profile, k)
            motifs.append(temp)
        if hamming_score(motifs) < hamming_score(bestMotifs):
            bestMotifs = motifs
        else:
           return bestMotifs



dna = [

'ACCCAGATTACGGGCGAGAGATTTCCAATGAATCATTTCAAGTCGCCACAAGGGTGCTCCTGAAAGGTGATGCACAGGGGTATACTGTCGATGCGAGGGCCTCGCGATCCGTATGGTCCGATAGACATTATTTAAATGCAGTAAGCCTGGAACCCAGATTACGGGC',
'GAGAGATTTCCAATGAATCATTTCAAGTCGCCACAAGGGTGCTCCTGAAAGGTGATGCACAGGGGTATACTGTCGATGCGAGGGCCTCGCGATCCGTATGGTCCGATAGACATTATTTAAATGCAGTAAGCCTGGAACCCAGCATCATTAAAAACGTATTACGGGC',
'CCGAGAAAATCTGGCACGGGACCGGAAGAGTGAGCTGGGAACTCTGCAGCAGCAAGTGCCAGAGATGCCAAGCGCTCATAGCACTCATCCTGACAGTCCCCTGACACATCTCTACGGGCTTAAGTAAGATCCACATAACAAAAAACGTACTGCTTCAACACCGCTC',
'AAATGCTCAATGCTCAATAAGTCCGCTACTTCGTGATGTCTAACCATCGAACTAACGTGCTCGGCATCTATCTGGGAACGAAAGTTCTTGAACGATAACTAGAAAAGCTCACACTAGGTCTGACTATGGCATAGGCGCAAACGTTTGTTAGAGGGTGTTAGCTATC',
'GGTTGTTCATAGGCTGAAACGTACGTACTCATTGACGAATTGCACTGCCATCTTAAGTCACCTCAGGCTCGCGCTTCCACACCCATGAAACTGAAAGGGGACTGGTTCTTCTCTTTTTAATCCGGCGGAGTAGTCTGGAATACATTCAGCCAGAAACATCGTCCCA',
'GGTACGGCGCAATGTGGGACGCCGCATGCCAGTGGTGAATGATTCCCTCTAATTGTTGTAGCAATGTTACACCGCTATGTAACTGGATAGGAAAGTCACGCTACAGGAAACGCTGTACAGATATCAGGCGTTAGGTAAAAACGGCACCTTCTACTTGCGGGGAAAA',
'CGCCGTTCCCACGTGTAAAAACGTTCTAATCACAGGTTAAGTAATTGTAGTTAGGAGTTGAAGACTGAGACTATCATGCATGAGATCTCATATAACTGTATAGCGGACCAGGCGCGTCTACAGCCTAGCTTTATGGGTGCCCCTGGCAACGACGCTCGAGCCCTGA',
'CTGCGTAAGGAACAAATCCGTTAACACATACTATAACAGTGAGCACGCCATGGCGGATAAACCTGCTCGAGGCGAATCCTGCATAGGTAGCGACGTCAACCCCCGTGATACCAGGGATATACGAAGCAGCGTTTGCTCGAGCTTGTGACTAACCGTGCCAAGATAC',
'CCCGGTCCAGCGTAATCGCTGGTGGCCAGCCTTCATAGAGAAACAGCACGTTTTCATAAGACAAACATAGGGACTATCATAGGTAAGCGCGTGTCAGCATATAATTTTCACAACGAATACCGACCGACGTTCGCACTCGACATACGGCGCCCGCGTCACACGCGGT',
'GTGCCAACTTGTACATAAGCAAGGTGTAGGATTATATTATTGACCCTCTCTGCAGAAGTTAATTCATAGAATAAAACGTTACGGAGGTGCGAGACAACAAAAAACTCAAAAAATGTGAGCTTCGCGCGCGGTTCCCCTTACTACTACCAGGTGACAAGAGGGATTT',
'ACTCCGCCATTTGGTTAAATAACTGCGGTCCCAGTCATTTAGAGGTATCGCAACATCTTGTGGTATGTTTAAGATCCCTACCTATCAAAGTCAGTATCTGTTCGGCTCGAGGAGCAGATAGGTAAAAACACGTTTCGTTAGAAAGCATTCATCCAGAGGGCAAAGG',
'CCAGGGGGCTCCACTTGGGTCGAGAGCATTTGGCTCATGAGAACCACAGACGGCATGACCCCGTCACGCTAAAGCGCTGTCAACCAGAAAGAAGACCGGGGTGTCATAGGTAAAGTTGTAGTGGATAGGTATTAAGGACTATATGAGTCGAGTCGTCCCCCGCTGC',
'GTTAAGGAGTGGGCCCTTAGGCCCTCCTTTGAAACTAAGCACTCGAGGACTGGTTGATGTCCTTTCTCTCAGTCAATGCACTTTCCAATAGGTCCTGCTTTGCTTGCTTTTCCTCTAGGAGGCTTTTTCGTCGTAACCCGTATTCATCCTTAAAAACGTAGGCAAG',
'GCAATGCGCCATACCGAAAAACGTGATGAGCACACAGAAGGGGAATCCGCGATCTGTACGGAGACCGTTAACCGTGCCCAAACACCGGTAGAGCGTGGGGTTCTGTACGCAGGAACACATAGATTTGCTGGACCTAGCGCTGCCTGTCCCCTGGACCTGCAATGGC',
'GGCTGTCCGAGACATTACAATTCGTGCTGCCCATCCGGTACGCGCCACGCCCGACGCTATAGATCAGCTCTTAAAATCCCCCCCCCGCGGGTCATTAGCGCTACGTTAGAAATCTATTTCCTTTTAACCATCAGGTAAAAACGTTACCCCGTATACCAAGTGCTCG',
'AGCTAAGAGATGCATAGGTCGGAACGTCACGATCCTGAGTCCCGTGTCATCTTGGCCGGGCCGGAGCTAGCGGATTAATCTACGTCACAATTTTAACGTATGAGGCGCCCACTCGCAGAGAATCAACATAGTACTGGGAATATAGATTCAACGTCTTGCAGAAGGA',
'CTAGGCCGCAGAGCCTCCATTTCATTGAGGATGGCTCGTTGCCTAAGCTGGACTCCTTGTAGGTACCGGGGTAAAAACGTCGTAATCGTGTCCCATCAAAGATACCATCCGCATCCTGTCTAGTGAGTCGTCGTAGATGGTTGTTTAATGAGGAAATTTCAACGGA',
'TGGTCGATGGGGACCACCATTAGCTATATCTTTCAAATCCTCTGCAGAAACACGGACGCGGCCGAGCCCGCCTCTTTAGCACGCATAGAAGAAAACGTACCCCTTCACGATGGAAGCTCGAGGGAGCTCACTTATTCCTATGCGGGGCGGGGTCATTGAATAATTC',
'CCAAGCGTAATCTCTGGTCGCTAGGGTTTGATGTAGGCTTACCTATCTTCCCATAAACTTTTATTGAGCCAAAGCCTATGGCTCCCGCTCGTCAAGTCTTTAGCGTGATCTGTAAAGTTGCATAGGTAAAAATTCTACGTCGGTAGGCTACATCTGCGGTCGTCAC',
'CGGGCGCTGTCCGGTCGTCACGACGGGTTACTGTATTTACGTCAGGAACTTCAGTCTTTCCGCAATTCCTAAGATAGTGTCTTGGAGCCGTAGTGTCAACTCTCACCGGCATGAACTCGGAAAAGCCCTCATAGGTAAAACTATGAGTCAGCGCTGGCCACACCAC',


]
k = 15
t = 20
num_of_trys = 10000
best_motifs = []

topMotifs = []
x = 0
while x < t:
    test = dna[x]
    topMotifs.insert(x, test[0:k])
    x += 1

current_dna = dna[0]

y = 0
go = 0
pseudocount = 1
beterMotifs = []
iii = 0
while iii < 1001:
    beterMotifs = randomized_motif_search(dna, k, t)
    if hamming_score(beterMotifs) < hamming_score(topMotifs):
        topMotifs = beterMotifs
        print(hamming_score(topMotifs))
    iii += 1
    print(iii)


for m in topMotifs:
    print(m)
print(hamming_score(topMotifs))





