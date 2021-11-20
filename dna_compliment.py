
def ReverseComplement(Pattern):   
    """
    Input:  A DNA string Pattern
    Output: The reverse complement of Pattern
    """
    Pattern = Reverse(Pattern) # reverse all letters in a string
    Pattern = Complement(Pattern) # complement each letter in a string
    return Pattern

def Reverse(Pattern):
    """
    First, reverse the dna pattern
    """
    return Pattern[::-1]


def Complement(Pattern):
    """
    Find the compliment the to the reversed strand
    """
    new_strand = ''
    d = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    for i in Pattern:
        for key in d:
            if i == key:
                new_word += d[key]
    return new_strand
  
  
def main():
    Pattern = 'AAAACCCGGT'
    compliment = ReverseComplement(Pattern)
    print(compliment)
  
  
if __name__ == "__main__":
    main()
