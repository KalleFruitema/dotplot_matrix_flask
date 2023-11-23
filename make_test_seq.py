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


def reverse():
    seq = random_seq(500)
    return seq, seq[::-1]


def inversie():
    seq1 = random_seq(500)
    seq2 = Seq(seq1[:200] + seq1[200:300][::-1] + seq1[300:])
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
        'reverse': reverse,
        'inverse': inversie, 
        'double_gap': double_gap, 
        'gap': gap}
    fasta_fmt = ''
    for naam, func in functies.items():
        seq1, seq2 = func()
        fasta_fmt += f">{naam}_1\n{seq1}\n>{naam}_2\n{seq2}\n"
    fasta_fmt = fasta_fmt.strip('\n')

    with open('test_sequences.fa', 'w') as file:
        file.write(fasta_fmt)


if __name__ == "__main__":
    main()
