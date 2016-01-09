#!/usr/bin/python
import argparse
import itertools

class HiddenMessages:
    def __init__(self):
        self.author = 'Jhalley de Castro'
        self.substitution = {
            'A': 'T',
            'T': 'A',
            'G': 'C',
            'C': 'G'
        }
        self.lex_order = {
            'A': 0,
            'C': 1,
            'G': 2,
            'T': 3
        }
        self.reverse_lex_order = {
            '0': 'A',
            '1': 'C',
            '2': 'G',
            '3': 'T'
        }

    def profile_most_probable_kmer(self, dna, k, pA, pC, pG, pT):
        profiles = {
            'A': pA,
            'C': pC,
            'G': pG,
            'T': pT
        }
        highest_prob = 0
        highest_prob_kmer = ''

        for i in xrange(len(dna) - k + 1):
            current_kmer = dna[i:i+k]
            current_kmer_prob = reduce(lambda a, b: a*b, [profiles[current_kmer[nuc_i]][nuc_i] for nuc_i in xrange(len(current_kmer))])

            if current_kmer_prob > highest_prob:
                highest_prob = current_kmer_prob
                highest_prob_kmer = current_kmer

        return highest_prob_kmer

    def median_string(self, dna, k):
        distance = 9999999999
        median = ''

        for pattern in itertools.product('ACTG', repeat=k):
            pattern = ''.join(pattern)
            d_pattern_dna = sum([min([self.hamming_distance(pattern, dna_i[j:j+k]) for j in range(len(dna_i) - k + 1)]) for dna_i in dna])

            if distance > d_pattern_dna:
                distance = d_pattern_dna
                median = pattern

        return median

    def motif_enumeration(self, dna_set, k, d):
        patterns = set()
        raw_patterns = set()
        for dna in dna_set:
            for raw_pattern in [dna[i:i+k] for i in range(len(dna) - k + 1)]:
                raw_patterns.add(raw_pattern)

        for pattern in raw_patterns:
            for pattern_prime in self.neighbors(pattern, d):
                appears = True
                for dna in dna_set:
                    if self.approx_pattern_count(pattern_prime, dna, d) == 0:
                        appears = False
                        break

                if appears:
                    patterns.add(pattern_prime)

        return patterns

    # From: http://stackoverflow.com/a/19941659
    def neighbors(self, word, hamming_distance, charset='ATCG'):
        for indices in itertools.combinations(range(len(word)), hamming_distance):
            for replacements in itertools.product(charset, repeat=hamming_distance):
                mutation = list(word)
                for index, replacement in zip(indices, replacements):
                    mutation[index] = replacement
                yield "".join(mutation)

    # This is the brute force approach
    def freq_words_w_mismatches_and_reverse_comp(self, text, k, d):
        freq = [0] * (4**k)

        # Calculate freq
        for i in xrange(4**k):
            pattern = self.number_to_pattern(i, k)
            freq[i] = self.approx_pattern_count(pattern, text, d) + self.approx_pattern_count(self.reverse_complement(pattern), text, d)

        # Return only max words
        maxFreq = max(freq)
        maxWords = []
        for i in xrange(4**k):
            if freq[i] == maxFreq:
                maxWords.append(self.number_to_pattern(i, k))

        return maxWords

    # This is the brute force approach
    def freq_words_w_mismatches(self, text, k , d):
        freq = [0] * (4**k)

        # Calculate freq
        for i in xrange(4**k):
            freq[i] = self.approx_pattern_count(self.number_to_pattern(i, k), text, d)

        # Return only max words
        maxFreq = max(freq)
        maxWords = []
        for i in xrange(4**k):
            if freq[i] == maxFreq:
                maxWords.append(self.number_to_pattern(i, k))

        return maxWords

    def approx_pattern_count(self, pattern, text, d):
        return len(self.approx_pattern_matching(pattern, text, d))

    def approx_pattern_matching(self, pattern, text, d):
        len_pattern = len(pattern)
        pos = []
        for i in xrange(len(text) - len_pattern + 1):
            if self.hamming_distance(text[i:i+len_pattern], pattern) <= d:
                pos.append(i)

        return pos

    def hamming_distance(self, p, q):
        hamming_distance = 0
        for i in xrange(len(p)):
            if p[i] != q[i]:
                hamming_distance += 1

        return hamming_distance

    def min_skew(self, genome):
        min_skew_value = 0
        min_skew_i = []
        skew_array = [0]

        for i in xrange(len(genome)):
            base = genome[i]
            if base == 'C':
                skew_array.append(skew_array[-1] - 1)
            elif base == 'G':
                skew_array.append(skew_array[-1] + 1)
            else:
                skew_array.append(skew_array[-1])

            if skew_array[-1] < min_skew_value:
                min_skew_value = skew_array[-1]
                min_skew_i = [i + 1]
            elif skew_array[-1] == min_skew_value:
                min_skew_i.append(i + 1)

        return min_skew_i 

    def skew_genome(self, genome):
        skew_array = [0]

        for base in genome:
            if base == 'C':
                skew_array.append(skew_array[-1] - 1)
            elif base == 'G':
                skew_array.append(skew_array[-1] + 1)
            else:
                skew_array.append(skew_array[-1])

        return skew_array

    def pattern_to_number2(self, pattern):
        if not pattern:
            return 0

        symbol = pattern[-1]
        prefix = pattern[:-1]
        return 4 * self.pattern_to_number2(prefix) + self.lex_order[symbol]

    def computing_frequencies(self, text, k):
        freq = [0] * (4**k)
        for i in xrange(len(text) - k + 1):
            pattern = text[i:i+k]
            j = self.pattern_to_number(pattern)
            #print '%s -> %s'%(pattern, j)
            freq[j] += 1

        return freq

    def number_to_pattern(self, index, k):
        # 11, 2 = GT
        def recurse(num):
            if num < 4:
                return str(num)
            else:
                return str(recurse(num/4)) + str(num%4)

        base4_in_dna = [self.reverse_lex_order[i] for i in recurse(index)]
        return ''.join(['A']*(k-len(base4_in_dna)) + base4_in_dna)
        

    def pattern_to_number(self, pattern):
        # This is essentially base 4
        # GT = 11
        number = 0
        for i in xrange(len(pattern)):
            number += 4**(len(pattern) - 1 - i) * self.lex_order[pattern[i]]

        return number

    def clump_finder2(self, genome, k, L, t):
        # Generate freq array 
        freq_array = {}

        for i in range(len(genome) - k + 1):
            if genome[i:i+k] not in freq_array:
                freq_array[genome[i:i+k]] = [i]
            else:
                freq_array[genome[i:i+k]].append(i)

        # Calculate Lt-clumps
        tWords = []
        for word in freq_array:
            if len(freq_array[word]) < t:
                continue
            else:
                for i in range(len(freq_array[word]) - t + 1):
                    if freq_array[word][i+t-1] + k - freq_array[word][i] <= L:
                        tWords.append(word)
                        break

        return tWords
        

    def clump_finder(self, genome, k, L, t):
        def frequent_words4(text, k, t):
            # Do everything in one pass
            wordCounts = {}
            tWords = {};
        
            for i in range(len(text) - k + 1):
                currWord = text[i:i+k]
                
                # Update freq count of word in currWords dict
                if currWord in wordCounts:
                    wordCounts[currWord] += 1
                else:
                    wordCounts[currWord] = 1
        
                # See whether this is a new tWord entry
                if wordCounts[currWord] >= t:
                    tWords[currWord] = True;
        
            return tWords.keys()

        kLt_patterns = []
        for i in xrange(len(genome) - L):
            window_genome = genome[i:i+L]
            new_patterns = frequent_words4(window_genome, k, t)
            kLt_patterns += new_patterns

        return set(kLt_patterns)

    def pattern_matching(self, pattern, genome):
        patternLength = len(pattern)
        positions = []
    
        for i in range(len(genome)):
            if genome[i:i+patternLength] == pattern:
                positions.append(i)
    
        return ' '.join([str(i) for i in positions])

    def reverse_complement(self, pattern):
        return ''.join([self.substitution[i] for i in pattern][::-1])

    def frequent_words(self, text, k):
        # Generate all words of len k from text
        words = {}
        for i in range(len(text) - k + 1):
            words[text[i:i+k]] = 0
    
        # Count freq of words
        # and add entry into wordsByFreq as well
        wordsByFreq = {}
        for word in words.keys():
            words[word] = self.pattern_count(text, word)
    
            if words[word] not in wordsByFreq:
                wordsByFreq[words[word]] = [word]
            else:
                wordsByFreq[words[word]].append(word)
    
        # Return maximum
        return ' '.join(wordsByFreq[max(wordsByFreq.keys())])
    
    def frequent_words2(self, text, k):
        # Do everything in one pass
        wordCounts = {}
    
        for i in range(len(text) - k + 1):
            currWord = text[i:i+k]
            if currWord in wordCounts:
                wordCounts[currWord] += 1
            else:
                wordCounts[currWord] = 1
    
        maxFreq = max(wordCounts.values())
        return ' '.join([i for i in wordCounts.keys() if wordCounts[i] == maxFreq])
    
    def frequent_words3(self, text, k):
        # Do everything in one pass
        wordCounts = {}
        currMax = 0;
        maxWords = [];
    
        for i in range(len(text) - k + 1):
            currWord = text[i:i+k]
    
            # Update freq count of word in currWords dict
            if currWord in wordCounts:
                wordCounts[currWord] += 1
            else:
                wordCounts[currWord] = 1
    
            # See whether this is a new maxWord entry
            if wordCounts[currWord] == currMax:
                maxWords.append(currWord)
            elif wordCounts[currWord] > currMax:
                currMax = wordCounts[currWord]
                maxWords = [currWord]
    
        return ' '.join(maxWords)

    def pattern_count(self, text, pattern):
        patternLength = len(pattern)
        count = 0
    
        for i in range(len(text)):
            if text[i:i+patternLength] == pattern:
                count += 1
    
        return count

if __name__ == "__main__":
    hm = HiddenMessages()

    parser = argparse.ArgumentParser(description='Coursera Hidden Messages - Jhalley')
    parser.add_argument('--pattern_count', help='Count the number of occurances of pattern in text')
    parser.add_argument('--freq_words', help='Find most frequent k-mers in text')
    parser.add_argument('--freq_words2', help='Find most frequent k-mers in text - faster')
    parser.add_argument('--freq_words3', help='Find most frequent k-mers in text - fastest!')
    parser.add_argument('--reverse_complement', help='Find the reverse complement of a DNA string')
    parser.add_argument('--pattern_matching', help='Find all occurrences of a pattern in a string')
    parser.add_argument('--clump_finder', help='Find patterns forming clumps in a string')
    parser.add_argument('--clump_finder2', help='Find patterns forming clumps in a string using freq array')
    parser.add_argument('--computing_frequencies', help='Generate frequency array')
    parser.add_argument('--pattern_to_number2', help='Pattern to number algorithm')
    parser.add_argument('--number_to_pattern', help='My number to pattern algorithm')
    parser.add_argument('--skew_genome', help='Get skew diagram of genome')
    parser.add_argument('--min_skew', help='Find a position in a genome where the skew diagram attains a minimum')
    parser.add_argument('--hamming_distance', help='Calculate hamming distance between two strings of equal length')
    parser.add_argument('--approx_pattern_matching', help='Find all approximate occurrences of a pattern in a string')
    parser.add_argument('--approx_pattern_count', help='Get count of all approximate occurrences of a pattern in a string')
    parser.add_argument('--freq_words_w_mismatches', help='Find the most frequent k-mers with mismatches in a string')
    parser.add_argument('--freq_words_w_mismatches_and_reverse_comp', help='Find the most frequent k-mers (with mismatches and reverse complements) in a string')
    parser.add_argument('--motif_enumeration', help='A brute force attempt at motif finding')
    parser.add_argument('--median_string', help='Find a median string')
    parser.add_argument('--profile_most_probable_kmer', help='Find a Profile-most probable k-mer in a string')
    args = parser.parse_args()

     
    # Parse from file
    if args.pattern_count:
        with open(args.pattern_count) as f:
            lines = f.readlines()
        print hm.pattern_count(lines[0].strip(), lines[1].strip())
    elif args.freq_words:
        with open(args.freq_words) as f:
            lines = f.readlines()
        print hm.frequent_words(lines[0].strip(), int(lines[1].strip()))
    elif args.freq_words2:
        with open(args.freq_words2) as f:
            lines = f.readlines()
        print hm.frequent_words2(lines[0].strip(), int(lines[1].strip()))
    elif args.freq_words3:
        with open(args.freq_words3) as f:
            lines = f.readlines()
        print hm.frequent_words3(lines[0].strip(), int(lines[1].strip()))
    elif args.reverse_complement:
        with open(args.reverse_complement) as f:
            lines = f.readlines()
        print hm.reverse_complement(lines[0].strip())
    elif args.pattern_matching:
        with open(args.pattern_matching) as f:
            lines = f.readlines()
        print hm.pattern_matching(lines[0].strip(), lines[1].strip())
    elif args.clump_finder:
        with open(args.clump_finder) as f:
            lines = f.readlines()
        klt = [int(i) for i in lines[1].strip().split()]
        print ' '.join(hm.clump_finder(lines[0].strip(), klt[0], klt[1], klt[2]))
    elif args.clump_finder2:
        with open(args.clump_finder2) as f:
            lines = f.readlines()
        klt = [int(i) for i in lines[1].strip().split()]
        print ' '.join(hm.clump_finder2(lines[0].strip(), klt[0], klt[1], klt[2]))
    elif args.computing_frequencies:
        with open(args.computing_frequencies) as f:
            lines = f.readlines()
        print ' '.join([str(i) for i in hm.computing_frequencies(lines[0].strip(), int(lines[1].strip()))])
    elif args.pattern_to_number2:
        with open(args.pattern_to_number2) as f:
            lines = f.readlines()
        print hm.pattern_to_number2(lines[0].strip())
    elif args.number_to_pattern:
        with open(args.number_to_pattern) as f:
            lines = f.readlines()
        print hm.number_to_pattern(int(lines[0].strip()), int(lines[1].strip()))
    elif args.skew_genome:
        with open(args.skew_genome) as f:
            lines = f.readlines()
        print hm.skew_genome(lines[0].strip())
    elif args.min_skew:
        with open(args.min_skew) as f:
            lines = f.readlines()
        print ' '.join([str(i) for i in hm.min_skew(lines[0].strip())])
    elif args.hamming_distance:
        with open(args.hamming_distance) as f:
            lines = f.readlines()
        print hm.hamming_distance(lines[0].strip(), lines[1].strip())
    elif args.approx_pattern_matching:
        with open(args.approx_pattern_matching) as f:
            lines = f.readlines()
        print ' '.join([str(i) for i in hm.approx_pattern_matching(lines[0].strip(), lines[1].strip(), int(lines[2].strip()))])
    elif args.approx_pattern_count:
        with open(args.approx_pattern_count) as f:
            lines = f.readlines()
        print hm.approx_pattern_count(lines[0].strip(), lines[1].strip(), int(lines[2].strip()))
    elif args.freq_words_w_mismatches:
        with open(args.freq_words_w_mismatches) as f:
            lines = f.readlines()
        print ' '.join(hm.freq_words_w_mismatches(lines[0].strip(), int(lines[1].strip().split()[0]), int(lines[1].strip().split()[1])))
    elif args.freq_words_w_mismatches_and_reverse_comp:
        with open(args.freq_words_w_mismatches_and_reverse_comp) as f:
            lines = f.readlines()
        print ' '.join(hm.freq_words_w_mismatches_and_reverse_comp(lines[0].strip(), int(lines[1].strip().split()[0]), int(lines[1].strip().split()[1])))
    elif args.motif_enumeration:
        with open(args.motif_enumeration) as f:
            lines = f.readlines()
        print ' '.join(sorted(hm.motif_enumeration([dna.strip() for dna in lines[1:]], int(lines[0].split()[0]), int(lines[0].split()[1]))))
    elif args.median_string:
        with open(args.median_string) as f:
            lines = f.readlines()
        print hm.median_string([dna.strip() for dna in lines[1:]], int(lines[0].strip()))
    elif args.profile_most_probable_kmer:
        with open(args.profile_most_probable_kmer) as f:
            lines = f.readlines()
        print hm.profile_most_probable_kmer(lines[0].strip(), int(lines[1].strip()), [float(i) for i in lines[2].strip().split()], [float(i) for i in lines[3].strip().split()], [float(i) for i in lines[4].strip().split()], [float(i) for i in lines[5].strip().split()])

    # Test calls
    #print hm.profile_most_probable_kmer('ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT', 5, [0.2, 0.2, 0.3, 0.2, 0.3], [0.4, 0.3, 0.1, 0.5, 0.1], [0.3, 0.3, 0.5, 0.2, 0.4], [0.1, 0.2, 0.1, 0.1, 0.2])
    #print hm.median_string(['AAATTGACGCAT', 'GACGACCACGTT', 'CGTCAGCGCCTG', 'GCTGAGCACCGG', 'AGTTCGGGACAG'], 3)
    #print sorted(hm.motif_enumeration(['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT'], 3, 1))
    # print set([i for i in hm.neighbors('TGCAT', 2)])
    # print len(set([i for i in hm.neighbors('TGCAT', 2)]))
    #print hm.freq_words_w_mismatches_and_reverse_comp('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1)
    # print hm.freq_words_w_mismatches('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1)
    # print len(hm.approx_pattern_matching('CCC', 'CATGCCATTCGCATTGTCCCAGTGA', 2))
    #print hm.approx_pattern_matching('ATTCTGGA', 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT', 3)
    # print hm.hamming_distance('CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA', 'CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG')
    #print hm.min_skew('TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT')
    # print hm.skew_genome('GCATACACTTCCCAGTAGGTACTG')
    #print hm.pattern_to_number2('ATGCAA')
    #print hm.computing_frequencies('ACGCGGCTCTGAAA', 2)
    #print hm.number_to_pattern(5437, 7)
    #print hm.number_to_pattern(5437, 8)
    #print hm.pattern_to_number('ATGCAA')
    #print hm.reverse_complement('CCAGATC')
    #print hm.pattern_matching('ATAT', 'GATATATGCATATACTT')
    #print hm.clump_finder('CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA', 5, 50, 4)
    #print hm.pattern_count('CGCGATACGTTACATACATGATAGACCGCGCGCGATCATATCGCGATTATC', 'CGCG')
    #print hm.frequent_words3('CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA', 3)




