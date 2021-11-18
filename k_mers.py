from collections import Counter
import more_itertools


def count_mers(Text, k):
    """
    Function takes in a string (Text) and mer length (k) and slides down the string finding repeating patterns
    returns a counter object containing a dictionary of pattern frequencies
    
    Notes:
    Uses 'sliding window' technique
    The window length is the value of 'k' and moves down the string one index at a time (step=1), creating k-length windows
    A Counter is a container that keeps track of how many times equivalent values are added
    
    """
    return Counter(("".join(mers) for mers in more_itertools.windowed(Text, k)))
  


def PatternCount(Text, Pattern):
    count = 0
    """
    Function takes in a string (Text) and a desired substring (Pattern) and finds the frequency of the that Pattern within the Text
    
    Notes: 
    We add the 1 as range is exclusive of the upper bound
    Thus, if ((10 - 3) + 1) = 8, we are moving over indexes 0 - 7
    Stopping at index 7 since we are checking a 'window' 3 indexes at a time (in this example)
    as 7, 8, 9 are the last three indexes of a len 10 string
    """
    for i in range(len(Text)-len(Pattern)+1):
        print(Text[i:i+len(Pattern)])
        # slicing indexes excludes upper bound
        # ex [0:3] == indexes 0, 1, 2
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count



def main():
    t = "ACAACTATGCATACTATCGGGAACTATCCT"
    # finds the frequency of a specific pattern
    c = PatternCount(t, "ACTAT")
    print(f"pattern count: {c}")

    # find all 5 mer length patterns
    d = count_mers(t, 5)
    print(f"patterns found: {d}")


if __name__ == "__main__":
    main()
