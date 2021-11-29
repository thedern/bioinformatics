def reverse_complement(pattern):
    """
    find the complement of a dna strand within a double-helix by first reversing the strand
    then finding its complement

    :param pattern: dna
    :type pattern: string
    :return: complementary dna string
    :rtype: string
    """

    pattern = pattern[::-1]
    new_strand = ''
    d = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    for i in pattern:
        for key, value in d.items():
            if i == key:
                new_strand += value
    return new_strand


def main():
    dna = 'AAAACCCGGT'
    compliment = reverse_complement(dna)
    print(compliment)


if __name__ == "__main__":
    main()
