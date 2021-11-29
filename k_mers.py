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

    count = 0
    # positions = []
    for i in range(len(text) - len(pattern) + 1):
        print(text[i:i + len(pattern)])
        # slicing indexes excludes upper bound
        # ex [0:3] == indexes 0, 1, 2
        if text[i:i + len(pattern)] == pattern:
            count += 1
            # positions.append(i)
    return count
    # return count, positions


def main():
    t = "ACAACTATGCATACTATCGGGAACTATCCT"
    # finds the frequency of a specific pattern
    c = pattern_counter("ACTAT", t)
    print(f"pattern count: {c}")

    # find all 5 mer length patterns
    d = frequency_map(t, 5)
    print(f"patterns found: {d}")


if __name__ == "__main__":
    main()
