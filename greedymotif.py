import sys
import re


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



dna = ["TAGAAGCTACGGCAAAGCTATCGCAATCAGACCAGTCTCCACTGTAACGATCTTTTCATAACGCCGGTCTTATTAGCTGAGCCGATTTAAAGTATACTGAGCATGCTGGTCTTAGGCTCGGCCACGATGCGTTTAGCCCCCATGAACGCTAACAAA",
"TTATGTAAGGTATGACACCCCATGAGACATTTCAAGACAGTGATTTCTATCGGTGTCGGTGTATACCGTATTCATTCCATAATCAGGAGAATGTATCTAGCGCCAACTGGTATTGCTTCAGGGCAACCCGCGCTTAGCTTCCATTTTTCCTATTGC",
"CCTTCCGGCCGAAGGTCCGGTAGGAAACATGGCCTGTCATACCCCATCTGTGTTCCTCTGACGCATGGTTCCCCTCTTTTGCCGCTTTACCGAGTGTCCTCGTTGGTTTTCTGTTGGCAGTGTGCCATCAGCAGGAAGACCAACTTAGATTCGTGT",
"CGATCCATAGTTTCGACCCTCAAGCCTAGGTTCCGTAGGGTGGTATTGTGTGGTGGGCTTGGATACGGACGCGTATTGCCGACCGTTATGTCGGCGCCCGGGATCAGCTGGACCGATAATGGGGCCGAATACGTCCACAGTTGCTAACTCCCCATG",
"ACCTCCTGTTTCACTTGCCTAACCTCAGCCCCCATTGGCGCAGATCACGCGAAATCAATGACAACGGATCATTAGTAAATCACGATCGTTCGCCCTGGACGAGGAGATATTACAGACTCATTACGAGCATTTCGTCATAATATATCAACCTGCACA",
"GTACCCGTGACAGCGCGAGGGAATGTCCGGCGAATCATTAGCGCGCGATCAGCCCCCATACATTAGAGCGGGAAACATGACCGTATCGTCCCTCTACCTAGGACCACTTTTTTTAGCCCCTTATCTTTGTCGAGCGACTTCCCGAGGGTTGACAAA",
"TGACGCCCCATCCTAGTGGAAATGTACAAGAAAGACTTCAGAATAAGGAGTGTCTATTAGGACGCGATGTAGAGTGTATGGGATCAAACTCGCATCCCGGCATATAGAGCGAAGTATGCACCGGGCCACTAGCGACGTCTCGAATCCGCCATGTCG",
"GGGCCAGTCGTCGACCCATTAGAGATGGACTTAAAATCAAACCCCATCCGCAGATGTTAGTTCTAGGAGAATATAGTAGACTGCAGCCAGGATCATCGCTCTCTTAGGGCTGGTATTGTGATTTGATGGGCTAAGTCCAGTTTACGAGGGCAAATA",
"CGCATAGGTCGAGATAGTCCGGCTGTTCTTGGAAATTTCTCGATCACAGAGTTAATTTGAAGTGTTTAGCTCGGATCTACTGTGTTAAGCCCCATGTTAACCCAGGTGCCCTTACATCTTCCGGTGCTGCCAGGGCGCGGTACAGCGTACGCCGTG",
"GATTGGATTCAGAGCCTCGACTAGGGTCTGCGCGACCCGGCCCCTCGATTACGCCCCATATTCCTTCCGGCCCCATGAGAGCCGAGGCGCGATTCTGAGATCGGACCCCCTACCATGTTCCGGCCAGCAGAAGAGGCAAATCACGCGAATCCCCGG",
"CATAACTCCCAGTAACTCCCCATTGGACCCACAGTGGTGGATAAGAGGACGCACTACGTCCCGGCACGGCGCTACCGTATGCCCAGACTGCCCACCGGTGTACCTCCGCAAGTACCATCCCGGTGTCAAATTGAGGTATCGGTATCAAATGCAGCC",
"ACATACAAGCGAATAAGAAGAGCATAATCCCCCATCTGTTAGGGTGCAAACTTTAGGCAAATGTAAAGGCGTGATTTAGGGACAATCTGGGACTTTGTTGTCTTGCGTGGAGCGCAACGACATCCTATGGAAAGACTAGCTACCTCCCTGGACGTC",
"TTGCTCTAACATTCAGACCCCATGTTACCGTAGAGACAGAGGGAAAACCCGGTAGTTACTTCTTCGCGGAAAAATGCGGATTCCCTTCCACCGCATGTGACCCTACCTCTCCCAGTCTTAGAACGACGTTTCATGGTTATCACATGGCAACGTGTG",
"TCACCTGATGACGGTTGCAATGACGTTTCCGCCACGTGATGCCCCATCGATATTAATATTTGCGCAGTACGTGAGTCAACCCCCGTGAAGAATCGTCATGATCTGAGTACTTGTTAGTGAAACTGCAGGAGAGTATAGCAGTGCTGTGCACCCTGA",
"GGTGGGAACTTAGAAAGCATACAGCCTACCGAGCTATACTTCAATTCCTCATTACGCGAGCAACGGGCGGGGTGGACCCGACCGGGGATCTTCTTGCTGTCAGACGTGTCACACCCCATTTGGCGAATCCCTCCACTCTCGACATCACAAATTCAA",
"AAGCCCTCGCAACTGAGTCGGGATGGTTAAGGACTGGGGAAGATTCTAAGCTTCTTGACTAGATCGGAAGGACGGTTCGAAGACTTGTATATGTCTAGGTTTTTAACCTAAGACCCCATAATGTTGGTTGTATGTTAATGGGACCAAACTGAGCGC",
"GTTTAATAGTTTAGATGCGTCCGTTAAGTCCCCATATGGCTAAGAGTACGGCAACCATTTACACCCGGATAAATGCGGGGAACAACCCAGGGAGAAAACGTCCACTAAACATGCTTTTAGCTTTATCATTAATACAAATTTTTCCTTAGCTACCGA",
"CGGTAGGTCGTGGAACCTGTATTTCAAAACGCCGACACCTGTAACATAATTTCACGCCCCTAATTGGGTCGGGGGGTTATTTTGTGACGGAGGTGAGATAGAGAAGCCGCCCGCTCTGGGTACAACCTCGATCCGTAGGATGAATCAGTCCCCATT",
"CACTCAATCTCGCGACCGTTGTCCACGTCCACTCGACCAGCTGAATGCTAAGGAACGTGCGACCATCCTACGGTTTTCTCATCATCACTCCCCATACCAAAAACCATGCGTCCCCTGTCCTGGATATTTTGACACTCTCAATAATCGACCGCTAGT",
"AAACCTATGAGCGGTGTCTGTTGGCCAATGGCCTCGGGTATCGAGTACTGCGGACCGAGCTGAAGCCCCATAAATACATTTTATAGAGTCCCAAGTCCAGGGGGACCTCCATACTACTCACGATCTCAAAGTGTCGCATATGCAAGGTTCGTCCTA",
"CTTCGGGGAACCCCGACAATGGCCTTATACCCCATCTTCGATCACCCGCTTCCATGATGATTGTGCTGTACTCAACCCGTCGTGGATTCTTTGCATGATTATTCCTGAGAACTGCTCGAGGTTCAGGTTGCCTGATCATTGAGACGGTAAGGAAAC",
"ATACCTTATTTCCCGAATTTACACGCATCTTCGGAATGCACAACATATACTTGAAGCCGGTAAGGCCCCATAAACACTTCTTGGTTCAGTCCAGCAGAGCAGACTCTCAGGACCGCCTAATCGTTAGGTCATTACTGGCTGGTCAAGAAGTGAGGG",
"GGCATTGAAAGGTTTAGTTCACGGGGTGGTCGCACTCAGCCGGTTCAAAGGACCAGTCTACATCATTCGCGGCTTACTTAAGGCCAAAGCCCAGAGTTTCATACAGAGCTGACGCAACTATAGATATGTACATCATCCCCCATAAAATGTCAGCGC",
"CGACTATGACAGCGCCCCCATTCCTCAAATGGGTATAAAAATTACTACATGCGTTAATTCGGTAGACGATTAAAGCGAAAAACGTGAGTCCCCATATACGCGCCATCAGGTCCGGTCGTTTTTTTCTCAAGTGGTCTAAGCGAACCAAGCCTTACC",
"CCGGAGAGACCGGTGCAGTGTATATATGCAGTTCGCCGATGATTCCAGACATCTTTAAGTTAGGAATCGTATTGCGTTTGTAAACAATTACCGCGGGGCTCCCTTTATGAGTAGGAGCAAGTTCTAGCATATTGAATCCCCATCTCTGTGCCATCC",
]

pseudocount = 1
k = 12
t = 25
i = 0
bestMotifs = []
x = 0
while x < t:
    test = dna[x]
    bestMotifs.insert(x, test[0:k])
    x += 1

current_dna = dna[0]
y = 0
while y < len(current_dna):
    kmer_motif = current_dna[y:(y + k)]
    if len(kmer_motif) < k:
        break
    count_profile = get_blank_profile(k, pseudocount)
    update_count_profile(kmer_motif, count_profile)
    probability_profile = calculate_prob_profile(count_profile, k)
    motifs = []
    motifs.append(kmer_motif)
    for l in range(1, t):
        temp_dna = dna[l]
        best_motife_of_dna = find_best_motif(temp_dna, count_profile, probability_profile, k)
        update_count_profile(best_motife_of_dna, count_profile)
        probability_profile = calculate_prob_profile(count_profile, k)
        motifs.append(best_motife_of_dna)
    if hamming_score(motifs) < hamming_score(bestMotifs):
        bestMotifs = motifs
    y += 1
for x in bestMotifs:
    print(x)









