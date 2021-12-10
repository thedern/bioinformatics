from collections import Counter
import more_itertools


def list_of_most_frequent_words(text, k):
    """
    Uses frequency_map_dictionary to return a list of most frequent words and their counts
    Gets the max value from the dict, and then captures all keys that have that value
    Note:  Also works with 'frequency_map' because counter object IS a dictionary (a dictionary subclass)

    :param text: dna sequence
    :type text: string
    :param k: length of mer
    :type k: int
    :return: list of most frequent mers
    :rtype: list
    """

    freq = frequency_map_dictionary(text, k)
    # freq = frequency_map(Text, k)
    m = max(freq.values())
    return [(key, m) for key in freq if freq[key] == m]


def frequency_map_dictionary(text, k):
    """
    This version of the frequency map produces the same net result as the code in 'frequency_map', but does not use intertools.
    The sliding window is produced vi index iteration.

    for example:  find all mers of length 5 in text

    :param text: dna sequence
    :type text: string
    :param k: length of mer
    :type k: int
    :return: dictionary of mers and the number of times they occur {'mer': count}
    :rtype: dict
    """

    freq = {}
    n = len(text)
    for i in range(n - k + 1):
        pattern = text[i:i + k]
        freq[pattern] = freq.get(pattern, 0) + 1
    return freq


def frequency_map(text, k):
    """
    Function takes in a string (Text) and mer length (k) and slides down the string finding repeating patterns
    returns a counter object containing a dictionary of pattern frequencies
    
    Notes:
    Uses 'sliding window' technique
    The window length is the value of 'k' and moves down the string one index at a time (step=1), creating k-length windows
    A Counter is a container that keeps track of how many times equivalent values are added
    
    Text = string
    K = int
    return = dict ('mers': counts)
    
    """
    return Counter(("".join(mers) for mers in more_itertools.windowed(text, k)))


def nucleotide_array(symbol, text):
    """
    Uses pattern_counter function to create a dict of index numbers and nucleotide counts {index: count}
    for a window of n-length starting at that index.

    for example: if the window length is 4, and the symbol is nucleotide 'a', the dict will return the number
    of 'a' nucleotides staring at index 0 for a window length of 4

    Loop example with a window of 4:
    first execution will be the number of 'a' staring at index 0, ending at index 3
    second execution will be the number of 'a' staring at index 1, ending at index 4
    third execution will be the number of 'a' staring at index 2, ending at index 5
    etc

    :param symbol: nucleotide (a, t, c, g)
    :type symbol: string
    :param text: dna sequence
    :type text: string
    :return: dictionary of indexes and counts
    :rtype: dict
    """

    n = len(text)
    extended_genome = text + text[0:n // 2]
    return {
        i: pattern_counter(symbol, extended_genome[i: i + (n // 2)])
        for i in range(n)
    }


def faster_nucleotide_array(symbol, text):
    """

    more efficient algorithm for a nucleotide array.

    :param symbol: nucleotide (a, t, c, g)
    :type symbol: string
    :param text: dna sequence
    :type text: string
    :return: dictionary of indexes and counts
    :rtype: dict
    """
    array = {}
    genome_len = len(text)
    extended_genome = text + text[0:genome_len // 2]

    """
    # look at the first half of Genome to compute first array value
    # Example: if 'A' is found 6 times in the half genome,  array key '0', (array[0])  = {0:6}
    """
    array[0] = pattern_counter(symbol, text[0:genome_len // 2])

    for i in range(1, genome_len):
        """
        # the current index value will always be based on the previous index value but can change as we do our analysis
        # we start our array with the number of symbols counted from position 0 for the length of the window
        # Example:  {0:6, 1:6}
        """
        array[i] = array[i - 1]

        """
        # get and analyze first nucleotide in window based on current index, which is beginning of window
        # if the first nucleotide is one we are looking for, we assume we loose it as we move one to the right by index
        # therefore, subtract
        # Example:  symbol 'A' window 3; [ATC]AGC => A[TCA]GC;  lost the first 'A' (minus -1)
        """
        beginning_nucleotide = extended_genome[i - 1]
        if beginning_nucleotide == symbol:
            array[i] = array[i] - 1

        """
        # get and analyze last nucleotide in window based on last index in window
        # as index moves to the right, it increases by 1, thus we need to -1 to maintain correct window size as we loop
        # if the nucleotide at the end of the window is one we are looking for, add it
        # Example:  symbol 'A' window 3; [ATC]AGC => A[TCA]GC;  gained an 'A' at end of window (add +1)
        """
        ending_nucleotide = extended_genome[i + (genome_len // 2) - 1]
        if ending_nucleotide == symbol:
            array[i] = array[i] + 1

        print(i, extended_genome[1:(genome_len // 2 + i)], symbol, ':', array[i])
    return array


def pattern_counter(pattern, text):
    """
    for any given pattern of nucleotides in the provided text, count how many times that pattern occurs
    optionally can return the position in the string where that pattern occurs

    :param pattern: nucleotide pattern
    :type pattern: text
    :param text: dna sequence
    :type text: string
    :return: count or count and position
    :rtype: int or tuple
    """

    return sum(
        text[i: i + len(pattern)] == pattern
        for i in range(len(text) - len(pattern) + 1)
    )


def skew_array(genome):
    """
    traversed from ori to ter in the 5' → 3' direction and are thus called forward half-strands
    traversed from ori to ter in the 3' → 5' direction and are thus called reverse half-strands
    Determine if on the forward or reverse half strand by determining the difference between G and C
    each forward half-strand has more guanine + 1
    each reverse half-strand has more cytosine - 1

    :param genome: dna
    :type genome: string
    :return: array of nucleotide counts
    :rtype: dict
    """
    a = {0: 0}
    for i in range(len(genome)):
        if genome[i] in ['A', 'T']:
            a[i + 1] = a[i]
        if genome[i] == 'G':
            a[i + 1] = a[i] + 1
        if genome[i] == 'C':
            a[i + 1] = a[i] - 1

    return list(a.values())


def minimum_skew(Genome):
    """
    locating ori: it should be found where the skew array attains a minimum.
    :param Genome: dna
    :type Genome: text
    :return: list positions where G is lowest, indicating ori
    :rtype: list
    """

    skew_list = skew_array(Genome)
    m = min(skew_list)
    return [i for i in range(len(skew_list)) if skew_list[i] == m]


def hamming_distance(p, q):
    # ham = 0
    # for x, y in zip(p, q):
    #     if x is not y:
    #         ham += 1
    # return ham
    return sum(x is not y for x, y in zip(p, q))


def approximate_pattern_matching(text, pattern, d):
    """
    Returns a list of positions (indexes) within the text where approximate pattern is found.
    Approximate pattern is any pattern less than or equal to the hamming distance of d
    Uses pattern_counter's algorithm

    :param pattern: nucleotide pattern
    :type pattern: dna
    :param text: string
    :type text: string
    :param d: max hamming distance
    :type d: int
    :return: list of indexes of possible pattern matches
    :rtype: list
    """

    return [
        i
        for i in range(len(text) - len(pattern) + 1)
        if hamming_distance(text[i: i + len(pattern)], pattern) <= d
    ]


def approximate_pattern_count(pattern, text, d):
    """
    Return the number of approximate pattern matches for pattern in text with hamming distance <= d
    Uses pattern_counter's algorithm

    :param pattern: nucleotide pattern
    :type pattern: dna
    :param text: string
    :type text: string
    :param d: max hamming distance
    :type d: int
    :return: number of possible pattern matches
    :rtype: int
    """
    return sum(
        hamming_distance(text[i: i + len(pattern)], pattern) <= d
        for i in range(len(text) - len(pattern) + 1)
    )


def main():
    t = "CATTCCAGTACTTCGATGATGGCGTGAAGA"
    # finds the frequency of a specific pattern
    # c = pattern_counter("ACTAT", t)
    # print(f"pattern count: {c}")

    a = 'TGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCAAACTCC'
    b = 'GAGCGATTAAGCGTGACAGCCCCAGGGAACCCACAAAACGTGATCGCAGTCCATCCGATCATACA'
    print(hamming_distance(a, b))


    # find all 5 mer length patterns
    # d = frequency_map(t, 5)
    # print(f"patterns found: {d}")

    # x = faster_nucleotide_array("A", t)
    # print(x)
    #
    x = minimum_skew(t)
    print(x)


if __name__ == "__main__":
    main()
