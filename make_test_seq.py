from Bio.Seq import Seq
from random import choice


def random_seq(len):
    return Seq(''.join(choice("ATCG") for _ in range(len)))


def complement():
    seq = random_seq(500)
    return seq, seq


def reverse_complement():
    seq = random_seq(500)
    return seq, seq.reverse_complement()


def inversie():
    seq1 = random_seq(500)
    seq2 = Seq(seq1[:200] + seq1.reverse_complement()[200:300] + seq1[300:])
    return seq1, seq2


def double_gap():
    seq1 = random_seq(500)
    seq2 = Seq(seq1[:200] + random_seq(100) + seq1[300:])
    return seq1, seq2


def gap():
    seq1 = random_seq(500)
    seq2 = Seq(seq1[:200] + seq1[300:])
    return seq1, seq2


def main():
    functies = {
        'comp': complement, 
        'rev_comp': reverse_complement, 
        'inverse': inversie, 
        'double_gap': double_gap, 
        'gap': gap}
    for naam, func in functies.items():
        seq1, seq2 = func()
        print(f"\n{naam}\n{seq1}\n{seq2}\n")


if __name__ == "__main__":
    main()
