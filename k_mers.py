from collections import Counter
import more_itertools


def list_of_most_frequent_words(Text, k):
    """
    Uses frequency_map_dictionary to return a list of most frequent words and their counts
    Gets the max value from the dict, and then captures all keys that have that value
    Note:  Also works with 'frequency_map' because counter object IS a dictionary (a dictionary subclass)
    
    Text = string
    K = int
    return = list of tuples [('mers', counts)]
    
    """
    words = []
    freq = frequency_map_dictionary(Text, k)
    # freq = frequency_map(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append((key, m))
    return words


def frequency_map_dictionary(Text, k):
    """
    This version of the frequency map produces the same net result as the code in 'frequency_map', but does not use intertools.
    The sliding window is produced vi index iteration.
    
    Text = string
    K = int
    return = dict ('mers': counts)
    
    """
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = freq.get(Pattern,0) + 1
    return freq


def frequency_map(Text, k):
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
    return Counter(("".join(mers) for mers in more_itertools.windowed(Text, k)))
  


def pattern_counter(Pattern, Text):
    count = 0
    """
    Function takes in a string (Text) and a desired substring (Pattern) and finds the frequency of the that Pattern within the Text
    
    Notes: 
    We add the 1 as range is exclusive of the upper bound
    Thus, if ((10 - 3) + 1) = 8, we are moving over indexes 0 - 7
    Stopping at index 7 since we are checking a 'window' 3 indexes at a time (in this example)
    as 7, 8, 9 are the last three indexes of a len 10 string
    
    Text = string
    Pattern =string
    count = int
    """
    positions = []
    for i in range(len(Text)-len(Pattern)+1):
        print(Text[i:i+len(Pattern)])
        # slicing indexes excludes upper bound
        # ex [0:3] == indexes 0, 1, 2
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
            positions.append(i)
    return count, positions



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
