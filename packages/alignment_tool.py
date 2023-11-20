"""
DISCLAIMER
Dit programma werkt voor de sequenties die ik voor BSEQAN heb
gekregen, maar ik kan niet garanderen dat dit ook werkt met 
andere sequenties. Sterker nog, ik weet wel zeker dat het niet 
altijd werkt. Maarja dit programma is toch voor de lol en ik heb
hier al te veel hoofdpijn van gekregen om het programma nog volledig
compatible te maken met alles wat je invoert. ¯\_(ツ)_/¯
"""


from pprint import pprint
from copy import deepcopy
from time import perf_counter


class Matrix:
    def __init__(self) -> None:
        pass

    def get_coords(self, column, row):
        column_val = self.matrix[0].index(column)
        for i, r in enumerate(self.matrix):
            if r[0] == row:
                row_val = i; break
        return column_val, row_val
    
    def get_val(self, column_i, row_i):
        return self.matrix[row_i][column_i]
    
    def get_score(self, column, row):
        return self.get_val(*self.get_coords(column.upper(), row.upper()))
    
    def read_csv(self, filepath):
        with open(filepath) as csv_file:
            matrix = csv_file.read().splitlines()
        for i, row in enumerate(matrix):
            matrix[i] = row.split(",")
        return matrix


class ScoreMatrix(Matrix):
    def __init__(self,  filepath: str) -> None:
        super().__init__()
        self.matrix = self.read_csv(filepath)
        self.gap_val = self.__find_gap_val()

    def __find_gap_val(self):
        for y in self.matrix:
            if y[0] == '-':
                return int(y[2])
            
    def get_gap_val(self):
        return self.gap_val
    

class AlignMatrix(Matrix):
    def __init__(self, seq1: str, seq2: str, scorematrix: ScoreMatrix, 
                 aligntype: str='N-W') -> None:
        self.seq1 = seq1
        self.seq2 = seq2
        self.aligntype = self.__check_aligntype(aligntype)
        self.scorematrix = self.__check_scorematrix(scorematrix)
        self.matrix = self.__make_matrix(seq1, seq2)
        self.__align()
        self.alignment, self.alignscore = self.__traceback()

    def __check_scorematrix(self, scorematrix):
        if type(scorematrix) == ScoreMatrix:
            return scorematrix
        else:
            raise AttributeError(f"{type(scorematrix)} is invalid for scorematrix, "
                                 "please use a \'ScoreMatrix\' class object")
    
    def __check_aligntype(self, aligntype):
        if aligntype in ['N-W', 'S-W']:
            return aligntype
        else:
            raise TypeError(f"\'{aligntype}\' is an invalid aligntype, choose either:\n"
                            "\'N-W\' for a global (full-sequence) alignment (default)\n"
                            "\'S-W\' for a local (partial-sequence) alignment")

    def __make_matrix(self, seq1, seq2):
        matrix = []
        for i in range(len(seq2) + 2):
            if i == 0:
                matrix.append([*f" -{seq1}"])
            elif i == 1:
                matrix.append((["-" if z == 0 else "" for z in range(len(seq1) + 2)]))
            else:
                matrix.append(([f"{seq2[i - 2]}" if z == 0 else "" for z in range(len(seq1) + 2)]))
        return matrix

    def __align(self):
        for ri, row in enumerate(self.matrix):
            for ci, col in enumerate(row):
                if ri == 0 or ci == 0:
                    continue
                if ri == 1 and ci == 1:
                    self.matrix[ri][ci] = 0
                elif ri == 1:
                    if self.aligntype == 'S-W':
                        self.matrix[ri][ci] = Score(0, None)
                    else:
                        self.matrix[ri][ci] = Score(self.matrix[ri][ci-1] + int(self.scorematrix.get_gap_val()), '→')
                elif ci == 1:
                    if self.aligntype == 'S-W':
                        self.matrix[ri][ci] = Score(0, None)
                    else:
                        self.matrix[ri][ci] = Score(self.matrix[ri-1][ci] + int(self.scorematrix.get_gap_val()), '↓')
                else:
                    up = self.matrix[ri-1][ci] + int(self.scorematrix.get_gap_val())
                    left = self.matrix[ri][ci-1] + int(self.scorematrix.get_gap_val())
                    diag_seq1 = self.matrix[0][ci]
                    diag_seq2 = self.matrix[ri][0]
                    diag = self.matrix[ri-1][ci-1] + int(self.scorematrix.get_score(diag_seq1, diag_seq2))
                    
                    best_direc = max([up, left, diag])
                    if self.aligntype == 'S-W':
                        best_direc = max([up, left, diag, 0])
                    if best_direc == up:
                        direction = '↓'
                    elif best_direc == left:
                        direction = '→'
                    elif best_direc == diag:
                        direction = '↘'
                    elif best_direc == 0 and self.aligntype == 'S-W':
                        direction = None
                    if direction != None:
                        if (up == left and direction in "↓→") or\
                            (up == diag and direction in "↓↘") or\
                                (left == diag and direction in "→↘"):
                            x = max(self.matrix[ri-1][ci], self.matrix[ri][ci-1], self.matrix[ri-1][ci-1])
                            if x == self.matrix[ri-1][ci]:
                                direction = '↓'
                            if x == self.matrix[ri][ci-1]:
                                direction = '→'
                            if x == self.matrix[ri-1][ci-1]:
                                direction = '↘'
                    self.matrix[ri][ci] = Score(best_direc, direction)

    def get_aligntable(self, origintable=False):
        table = ""
        for ri, row in enumerate(self.matrix):
            for ci, col in enumerate(row):
                try:
                    if origintable:
                        if col.get_origin() == None:
                            table += (f"{'0':>4}")
                        else:
                            table += (f"{col.get_origin():>4}")
                    else:
                        table += (f"{col.__str__():>4}")
                except:
                    table += (f"{col:>4}")
                if ci == 0:
                    table += (f"  |")
            table += '\n'
            if ri == 0:
                table += f"{(len(row) + 1)*'----'}\n"
        return table
    
    def __traceback(self):
        row = col = -1
        if self.aligntype == 'S-W':
            max_list = []
            for irow, m_row in enumerate(self.matrix):
                for icol, col_val in enumerate(m_row):
                    if type(col_val) in [int, Score]:
                        max_list.append((col_val, irow, icol))
            start = max(max_list, key=lambda x: x[0])
            row, col = start[1], start[2]

        alignment = [[self.matrix[0][col]],[],[self.matrix[row][0]]]
        totalscore = 0
        pos = self.matrix[row][col]
        totalscore += pos.get_number()
        while True:
            try:
                match pos.get_origin():
                    case '↖':
                        row -= 1
                        col -= 1
                        pos = self.matrix[row][col]
                        alignment[0].append(self.matrix[0][col])
                        alignment[2].append(self.matrix[row][0])
                        totalscore += pos.get_number()
                        if pos == 0:
                            break
                    case '↑':
                        row -= 1
                        pos = self.matrix[row][col]
                        alignment[0].append('-')
                        alignment[2].append(self.matrix[row][0])
                        totalscore += pos.get_number()
                        if pos == 0:
                            break
                    case '←':
                        col -= 1
                        pos = self.matrix[row][col]
                        alignment[0].append(self.matrix[0][col])
                        alignment[2].append('-')
                        totalscore += pos.get_number()
                        if pos == 0:
                            break
                    case None:
                        break
            except AttributeError:
                break
        alignment[0] = alignment[0][:-1]
        alignment[2] = alignment[2][:-1]
        for i, am in enumerate(alignment.copy()):
            if i == 1:
                continue
            for bi, base in enumerate(am.copy()):
                if base == '-':
                    alignment[i][bi], alignment[i][bi-1] = alignment[i][bi-1], base
        alignment[0] = alignment[0][::-1]
        alignment[2] = alignment[2][::-1]
        for i, base in enumerate(alignment[0]):
            if base == alignment[2][i]:
                alignment[1].append('|')
            elif base == '-' or alignment[2][i] == '-':
                alignment[1].append(' ')
            else:
                alignment[1].append('*')
        return alignment, totalscore
    
    def get_alignment(self):
        return self.alignment
    
    def get_real_alignscore(self, gap_open: int=-10, gap_ext: int=-1,
                            scorematrix: ScoreMatrix=None):
        if scorematrix == None:
            scorematrix = self.scorematrix
        real_score = 0
        gap = False
        for i, char in enumerate(self.alignment[1]):
            if char in "|*":
                gap = False
                real_score += int(scorematrix.get_score(self.alignment[0][i],
                                                        self.alignment[2][i]))
            else:
                if gap == False:
                    real_score += gap_open
                    gap = True
                else:
                    real_score += gap_ext
        return real_score
    
    def get_alignscore(self):
        return self.alignscore
    
    def __str__(self) -> str:
        text = ''
        for i in self.alignment:
            text += f'{" ".join(i)}\n'
        return text
    

class Score:
    def __init__(self, number, direction) -> None:
        self.number = number
        self.direction = direction
        self.origin = self.__make_origin()

    def __make_origin(self):
        match self.direction:
            case '↓':
                return '↑'
            case '→':
                return '←'
            case '↘':
                return '↖'
            
    def get_number(self):
        return self.number
            
    def get_origin(self):
        return self.origin

    def __eq__(self, __value: object) -> bool:
        return __value == self.number
    
    def __add__(self, other):
        return self.number + other
    
    def __sub__(self, other):
        return self.number - other
    
    def __lt__(self, other):
        try:
            return self.number < other.number
        except AttributeError:
            return self.number < other
    
    def __gt__(self, other):
        try:
            return self.number > other.number
        except AttributeError:
            return self.number > other
    
    def __repr__(self) -> str:
        return str(self.number)


class DistanceMatrix(Matrix):
    def __init__(self, sequences: str=None, scorematrix: ScoreMatrix=None, matrix_filepath: str=None) -> None:
        super().__init__()
        self.scorematrix = scorematrix
        if sequences == None == matrix_filepath:
            raise ValueError("Neither \'sequences\' and \'matrix_filepath\' were given, while exactly one is required.")
        elif None not in [sequences, matrix_filepath]:
            raise ValueError("Both \'sequences\' and \'matrix_filepath\' were given while exactly one is allowed.")
        if sequences != None: 
            sequences = self.__make_sequence_dict(sequences)
            alignment_dict = self.__get_alignments(sequences)
            self.distmatrix = self.__make_matrix(sequences, alignment_dict)
        else:
            self.distmatrix = self.read_csv(matrix_filepath)
            for i, y in enumerate(self.distmatrix):
                for x, z in enumerate(y):
                    if z.isnumeric():
                        try:
                            self.sequences[i][x] = int(x)
                        except Exception:
                            continue

    def __read_sequences(self, filepath):
        with open(filepath) as file:
            inhoud = file.read().strip().split('\n')
        seq_dict = {}
        for i, line in enumerate(inhoud):
            if line[0] == '>':
                seq_dict.update({
                    line: inhoud[i+1]
                })
        return seq_dict

    def __get_alignments(self, sequences):
        alignment_dict = {}
        for header, seq in sequences.items():
            alignment_dict.update({
                header : {}
            })
            for header2, seq2 in sequences.items():
                alignment = AlignMatrix(seq, seq2, self.scorematrix, aligntype="S-W")
                score = alignment.get_real_alignscore()
                alignment_dict[header].update({
                    header2: score
                })
        return alignment_dict

    def __make_sequence_dict(self, sequences):
        try:
            sequences = self.__read_sequences(sequences)
        except FileNotFoundError:
            if type(sequences) == list:
                seq_copy = sequences.copy()
                sequences = {}
                for i, seq in enumerate(seq_copy):
                    sequences.update({
                        f">seq_{i+1:>03}": seq
                    })
            elif (seq_type := type(sequences)) != dict:
                raise ValueError(f"{seq_type} used for sequences, possible types are:\n"\
                                 "string: ONLY if inputting file location for .txt .fa or .fna file\n"\
                                 "list: list of (string) sequences (headers will be generated)\n"\
                                 "dict: headers (string) paired with sequences (string)")
        if (seq_len := len(sequences)) < 3:
            raise ValueError(f"Only {seq_len} sequence(s) given, minumum is 3.")
        return sequences
    
    def __make_matrix(self, sequences, alignment_dict):
        dist_matrix = []
        for i, header in enumerate(sequences):
            if i == 0:
                dist_matrix.append([header for header in sequences])
                dist_matrix[0].insert(0, '')
            dist_matrix.append([header, *['' for _ in sequences]])
        for ri, row in enumerate(dist_matrix):
            if ri == 0:
                continue
            for ci, col in enumerate(row):
                if ci == 0:
                    continue
                head1, head2 = dist_matrix[0][ci], dist_matrix[ri][0]
                self_align1, self_align2 = alignment_dict[head1][head1], alignment_dict[head2][head2]
                cross_align = alignment_dict[head1][head2]
                dist_matrix[ri][ci] = (1 - cross_align / self_align1) * (1 - cross_align / self_align2)
        return dist_matrix
    
    def find_lowest_score(self, matrix: list[list]=None):
        lowest = []
        if matrix == None:
            matrix = self.distmatrix
        for ri, row in enumerate(matrix):
            if ri == 0:
                continue
            lowest.append(sorted(row, key=lambda x: x if type(x) != str else 1)[1])
        lowest = min(lowest)
        breaker = False
        for ri, row in enumerate(matrix):
            if ri == 0:
                continue
            for ci, val in enumerate(row):
                if val == lowest:
                    lowest_seqs = (matrix[0][ci], matrix[ri][0])
                    breaker = True
                    break
            if breaker:
                break
        return lowest, *lowest_seqs
    
    def __make_new_matrix(self, matrix, tree_vals=None, prints: bool=False):
        t_matrix = deepcopy(matrix)
        if tree_vals is None:
            tree_vals = []
        lowest, seq1, seq2 = self.find_lowest_score(t_matrix)
        curr_pair = OTUpair(lowest, seq1, seq2=seq2)
        tree_vals.append(curr_pair)
        n_matrix = []
        top_added = False
        left_added = False
        skip_line = False
        for ri, row in enumerate(t_matrix):
            line = []
            for ci, val in enumerate(row):
                if ri == 0:
                    if val in (seq1, seq2):
                        if not top_added:
                            line.append(curr_pair.__repr__())
                            top_added = True
                    else:
                        line.append(val)
                elif ci == 0:
                    if val in (seq1, seq2):
                        if not left_added:
                            line.append(curr_pair.__repr__())
                            left_added = True
                        else:
                            skip_line = True
                    else:
                        line.append(val)
                else:
                    line.extend(["" for _ in n_matrix[0][:-1]])
                    break
                    
            if not skip_line:
                n_matrix.append(line)
            skip_line = False
        if prints:
            print(self.format_matrix(n_matrix))
        for ri, row in enumerate(n_matrix):
            if ri == 0:
                continue
            for ci, val in enumerate(row):
                if ci == 0:
                    continue
                if n_matrix[ri][0] == n_matrix[0][ci]:
                    n_matrix[ri][ci] = 0
                elif "," in (comb := n_matrix[0][ci]) and comb not in t_matrix[0]:
                    comb_pair = [pair for pair in tree_vals if pair.comb_seq == comb][0]
                    left_idx = [i for i, y in enumerate(t_matrix) if y[0] == n_matrix[ri][0]][0]
                    dist1 = t_matrix[t_matrix[0].index(comb_pair.seq1)][left_idx]
                    dist2 = t_matrix[t_matrix[0].index(comb_pair.seq2)][left_idx]
                    comb_dist = (dist1 + dist2) / 2
                    n_matrix[ri][ci] = comb_dist 
                elif "," in (comb := n_matrix[ri][0]) and comb not in t_matrix[0]:
                    comb_pair = [pair for pair in tree_vals if pair.comb_seq == comb][0]
                    left_idx1 = [i for i, y in enumerate(t_matrix) if y[0] == comb_pair.seq1][0]
                    left_idx2 = [i for i, y in enumerate(t_matrix) if y[0] == comb_pair.seq2][0]
                    dist1 = t_matrix[t_matrix[0].index(n_matrix[0][ci])][left_idx1]
                    dist2 = t_matrix[t_matrix[0].index(n_matrix[0][ci])][left_idx2]
                    comb_dist = (dist1 + dist2) / 2
                    n_matrix[ri][ci] = comb_dist 
                else:
                    left_idx = [i for i, y in enumerate(t_matrix) if y[0] == n_matrix[ri][0]][0]
                    comb_dist = t_matrix[t_matrix[0].index(n_matrix[0][ci])][left_idx]
                    n_matrix[ri][ci] = comb_dist
                if prints:
                    print(self.format_matrix(n_matrix))

        return n_matrix, tree_vals

    def __newick_format(self, matrix, tree_vals) -> str:
        remaining_vals = [(i, val) for i, val in enumerate(matrix[0]) if val != '']
        for r_val in remaining_vals:
            if r_val[1] not in (obj.comb_seq for obj in tree_vals):
                left_idx = [i for i, y in enumerate(matrix) if y[0] not in ("", r_val[1])][0]
                dist = matrix[r_val[0]][left_idx]
                tree_vals.append(OTUpair(dist, r_val[1], outgroup=True))
        tree_str = ''
        check_tree_vals = [i.__repr__() for i in tree_vals]

        # for tree_val in tree_vals:
        #     for tree_val2 in tree_vals:
        #         inside_list = []
        #         if tree_val.comb_seq == tree_val2.comb_seq:
        #             continue
        #         if tree_val2.comb_seq in tree_val.comb_seq:
        #             inside_list.append(tree_val2)            
        #     tree_notation = ''
        #     inside_list.sort(key= lambda x: tree_val)
        #     for val in inside_list:
        #         tree_notation += f"({val.seq1}:{val.distance/2},{val.seq2}:{val.distance/2}):{tree_val.distance/2 - val.distance/2}"
        #     tree_str += tree_notation

        # for i in tree_vals:
        #     if i.outgroup:
        #         tree_notation

        #         break
                
        for tree_val in tree_vals:
            check = [i for i in tree_vals if i.__repr__() != tree_val.__repr__() and tree_val.__repr__() in i.__repr__()]
            if check != []:
                check = check[0]
                if check.__repr__() in check_tree_vals:
                    tree_str += check.tree_notation
                    check_tree_vals.pop(check_tree_vals.index(check.__repr__()))
                    
                if check.__repr__().index(tree_val.__repr__()) == 0:
                    tree_str = tree_str.replace("1?", f"{round(check.distance/2 - tree_val.distance/2, 2)}")
                else:
                    tree_str = tree_str.replace("2?", f"{round(check.distance/2 - tree_val.distance/2, 2)}")
                tree_str = tree_str.replace(tree_val.__repr__(), tree_val.tree_notation)
                check_tree_vals.pop(check_tree_vals.index(tree_val.__repr__()))

        if len(check_tree_vals) != 0:
            missed = [i for i in tree_vals if i.__repr__() in check_tree_vals][0]
            x, y = tree_str[::-1].split(":", 1)
            x = round(missed.distance/2 - float(x[::-1].strip()), 2)
            y = y[::-1]
            tree_str = f"{y}:{x}"
            tree_str += f",{missed.tree_notation}"
        tree_str = f"({tree_str});"

        return tree_str

    def make_tree(self, outputfile: str="tree_1.ph", prints: bool=False):
        matrix = self.distmatrix
        tree_vals = []
        while len(matrix) > 3:
            matrix, tree_vals = self.__make_new_matrix(matrix, tree_vals)
            if prints:
                print(f"\n{self.format_matrix(matrix)}\n\n")
        newick_tree = self.__newick_format(matrix, tree_vals)
        if prints:
            print(f"\nNewick format:\n{newick_tree}\n")
        if "." not in outputfile:
            outputfile += ".ph"
        with open(outputfile, "w+") as file:
            file.write(newick_tree)

    def get_distance_matrix(self):
        return self.distmatrix
    
    def format_matrix(self, matrix) -> str:
        fmt_matrix = deepcopy(matrix)
        fmt_matrix.insert(1, ["-"*11 for _ in range(len(fmt_matrix[0]) + 1)])
        for i, row in enumerate(fmt_matrix):
            if i == 1:
                continue
            fmt_matrix[i].insert(1, " | ")
        for ri, row in enumerate(fmt_matrix):
            for ci, val in enumerate(row):
                if type(val) == float:
                    fmt_matrix[ri][ci] = round(val, 3)
        str_matrix = "\n".join(["".join([f"{val:>11}" for val in line]) for line in fmt_matrix])
        return str_matrix
    
    def __str__(self) -> str:
        return self.format_matrix(self.distmatrix)


class OTUpair:
    def __init__(self, distance, seq1, seq2=None, outgroup=False) -> None:
        self.distance = distance
        self.seq1 = seq1
        self.outgroup = outgroup
        if not outgroup:
            self.seq2 = seq2
            self.comb_seq = f"{seq1},{seq2}"
            if self.comb_seq.count(",") == 1:
                self.tree_notation = f"({seq1}:{round(distance/2, 2)},{seq2}:{round(distance/2, 2)})"
            else:
                self.tree_notation = f"({seq1}:1?,{seq2}:2?):{round(self.distance/2, 2)}"
        else:
            self.comb_seq = f"{seq1}"
            self.tree_notation = f"{seq1}:{round(distance/2, 2)}"

    def __str__(self) -> str:
        return self.comb_seq

    def __repr__(self) -> str:
        return self.comb_seq


def main():
    matrix_path = "kimura_matrix.csv"
    seq1 = "ACCGAATTTTGCT"
    seq2 = "ACTCGAATCCT"
    # seq1 = "ATCCT"
    # seq2 = "ATTGACCA"
    nucl_score_matrix = ScoreMatrix(matrix_path)
    alignment = AlignMatrix(seq1, seq2, nucl_score_matrix)
    print(alignment.get_aligntable())
    print(alignment.get_aligntable(origintable=True))
    print(alignment)
    print(alignment.get_alignscore())

    seq3 = 'GKLEWTCAGHN'
    seq4 = 'ADGWERTGCA'
    aa_score_matrix = ScoreMatrix('PAM300.csv')
    aa_alignment = AlignMatrix(seq3, seq4, aa_score_matrix)
    print(aa_alignment.get_aligntable())
    print(aa_alignment.get_aligntable(origintable=True))
    print(aa_alignment)
    print(aa_alignment.get_alignscore())

    aa_alignment2 = AlignMatrix(seq3, seq4, aa_score_matrix,
                                aligntype='S-W')
    print(aa_alignment2.get_aligntable())
    print(aa_alignment2.get_aligntable(origintable=True))
    print(aa_alignment2)
    print(aa_alignment2.get_alignscore())

    seq5 = 'CGGACGGATGACC'
    seq6 = 'AACTAGACGCTCTTAAGT'
    alignment2 = AlignMatrix(seq5, seq6, nucl_score_matrix,
                             aligntype='S-W')
    print(alignment2.get_aligntable())
    print(alignment2.get_aligntable(origintable=True))
    print(alignment2)
    print(alignment2.get_alignscore())
    
    print(alignment2.get_real_alignscore())

    dm = DistanceMatrix(sequences="tree_seqs.fa", scorematrix=nucl_score_matrix)
    print(dm)


def main2():
    matrix_path = "kimura_matrix.csv"
    nucl_score_matrix = ScoreMatrix(matrix_path)

    dm = DistanceMatrix(sequences="tree_seqs_test.fa", scorematrix=nucl_score_matrix)
    print(dm)


def main_opdracht():
    blosum80 = ScoreMatrix('BLOSUM80.csv')

    t1 = perf_counter()
    dm2 = DistanceMatrix(sequences="opdracht_seqs.fa", scorematrix=blosum80)
    t2 = perf_counter()
    print(dm2)
    print(f"\n\nTook {t2-t1} seconds.")

    dm2.make_tree(outputfile="bseqan_phylo_tree.ph", prints=True)

    t3 = perf_counter()
    dm3 = DistanceMatrix(sequences="opdracht_seqs_namen.fa", scorematrix=blosum80)
    t4 = perf_counter()
    
    dm3.make_tree(outputfile="bseqan_phylo_tree_namen.ph")


if __name__ == "__main__":
    main_opdracht()
