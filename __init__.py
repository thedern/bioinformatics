
from k_mers_class import KMERS


def main():
    m = KMERS()
    t = "ACAACTATGCATACTATCGGGAACTATCCT"
    # finds the frequency of a specific pattern
    c = m.pattern_counter("ACTAT", t)
    print(f"pattern count: {c}")


if __name__ == "__main__":
    main()
