from collections import Counter
import more_itertools


class KMERS:
    def __init__(self, dna_sequence, nucleotide, mer_length=None):
        self.k = mer_length
        self.text = dna_sequence
        self.pattern = nucleotide

    @property
    def extended_genome(self):
        # automatically set extended_genome if len of pattern is longer than a single nucleotide
        n = len(self.text)
        if n > 1:
            return self.text + self.text[0:n // 2]

    def frequency_map(self):
        """
        Function takes in a string (Text) and mer length (k) and slides down the string finding repeating patterns
        returns a counter object containing a dictionary of pattern frequencies

        Notes:
        Uses 'sliding window' technique
        The window length is the value of 'k' and moves down the string one index at a time (step=1), creating k-length windows
        A Counter is a container that keeps track of how many times equivalent values are added

        """
        if self.k:
            return Counter(("".join(mers) for mers in more_itertools.windowed(self.text, self.k)))
        else:
            return "please set mer length"

    def list_of_most_frequent_words(self):
        """
        Uses frequency_map_dictionary to return a list of most frequent words and their counts
        Gets the max value from the dict, and then captures all keys that have that value
        Note:  Also works with 'frequency_map' because counter object IS a dictionary (a dictionary subclass)

        :return: list of most frequent mers
        :rtype: list
        """

        freq = self.frequency_map()
        # freq = frequency_map(Text, k)
        m = max(freq.values())
        return [(key, m) for key in freq if freq[key] == m]

    def pattern_counter(self, text=None):
        """
        for any given pattern of nucleotides in the provided text, count how many times that pattern occurs
        optionally can return the position in the string where that pattern occurs

        :param text: self.text or extended_genome
        :type: string
        :return: count or count and position
        :rtype: int or tuple
        """

        if text is None:
            text = self.text

        return sum(
            text[i: i + len(self.pattern)] == self.pattern
            for i in range(len(text) - len(self.pattern) + 1)
        )

    def nucleotide_array(self):
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

        :return: dictionary of indexes and counts
        :rtype: dict
        """
        n = len(self.text)

        return {
            i: self.pattern_counter(self.extended_genome[i: i + (n // 2)])
            for i in range(n)
        }
