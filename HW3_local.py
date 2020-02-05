import sys


def local_alignment(seq_1, seq_2, score_tuple):
    '''
        Finds the optimal local alignment between two sequences, using
        a scoring tuple structured as follows:

        score_tuple: [match, mismatch, gap-extend, *gap-open]
        *ignored
    '''
    D = []
    alignment_table = [] #"diag", "vert", "hori", "zero"

    seq_1_padded = "-" + seq_1
    seq_2_padded = "-" + seq_2

    for x in range(len(seq_1)+1):
        D.append([False]*(len(seq_2)+1))
        alignment_table.append(["zero"]*(len(seq_2)+1))

    for i in range(len(seq_1)+1):
        D[i][0] = 0

    for j in range(len(seq_2)+1):
        D[0][j] = 0

    for i in range(1,len(seq_1)+1):
        for j in range(1,len(seq_2)+1):
            recur_helper(D, i, j, seq_1_padded[i], seq_2_padded[j], score_tuple, alignment_table)

    best_pos = (-1, -1)
    best_score = -99
    for i in range(0,len(seq_1)+1):
        for j in range(0,len(seq_2)+1):
            if D[i][j] > best_score:
                best_score = D[i][j]
                best_pos = (i,j)

    best_align = TraceBack(best_pos, alignment_table, seq_1_padded, seq_2_padded)

    return best_score, best_align

def recur_helper(D, i, j, cur_char_1, cur_char_2, score_tuple, alignment_table):
    '''
        Performs the table filling for a single cell in the D table
        as well as the alignment table.
    '''
    match_or_mismatch = 0 #assume match
    if cur_char_1 != cur_char_2:
        match_or_mismatch = 1 #it's a mismatch

    diag_score = (D[i-1][j-1] + score_tuple[match_or_mismatch])
    vert_score = (D[i-1][j] + score_tuple[2])
    hori_score = (D[i][j-1] + score_tuple[2])

    if diag_score >= vert_score and diag_score >= hori_score and diag_score >= 0:
        alignment_table[i][j] = "diag"
        D[i][j] = diag_score
    elif vert_score >= diag_score and vert_score >= hori_score and vert_score >= 0:
        alignment_table[i][j] = "vert"
        D[i][j] = vert_score
    elif hori_score >= diag_score and hori_score >= vert_score and hori_score >= 0:
        alignment_table[i][j] = "hori"
        D[i][j] = hori_score
    else:
        alignment_table[i][j] = "zero"
        D[i][j] = 0

def TraceBack(best_pos, alignment_table, seq_1, seq_2):
    '''
        Traces back the path given in the alignment table to
        create a sequence alignment
    '''
    align1 = ""
    align2 = ""

    cur_pos = best_pos
    while alignment_table[cur_pos[0]][cur_pos[1]] != "zero":
        if alignment_table[cur_pos[0]][cur_pos[1]] == "vert":
            align1 = seq_1[cur_pos[0]] + align1
            align2 = "-" + align2
            cur_pos = (cur_pos[0]-1, cur_pos[1])
        elif alignment_table[cur_pos[0]][cur_pos[1]] == "hori":
            align1 = "-" + align1
            align2 = seq_2[cur_pos[1]] + align2
            cur_pos = (cur_pos[0], cur_pos[1]-1)
        else: #diag
            align1 = seq_1[cur_pos[0]] + align1
            align2 = seq_2[cur_pos[1]] + align2
            cur_pos = (cur_pos[0]-1, cur_pos[1]-1)

    return align1, align2


def FASTAparse(filename):
    '''
        Reads a FASTA file and turns it into a single string.
    '''
    file_lines = ""
    file = open(filename, "r")
    file.readline()

    for line in file:
        nextline = line.strip()
        file_lines = file_lines + nextline

    file_lines.lower()
    file.close()
    return file_lines

def scoring_parse(filename):
    '''
        Reads a scoring file and turns it into an integer array.
    '''
    with open(filename, "r") as file:
        file.readline()
        score_tuple = file.readline().split()
        return [int(i) for i in score_tuple]

def printNice(best_score, best_align):
    '''
        Takes an optimal score and a best alignment pair,
        and prints them in a readable fashion.
    '''
    print("Alignment score:", best_score)
    print("Optimal Local Alignment:")
    print(best_align[0])
    print(best_align[1])

def main():
    FASTA_file_1 = sys.argv[1]
    FASTA_file_2 = sys.argv[2]
    scoring_file = sys.argv[3]

    seq1_lines = FASTAparse(FASTA_file_1)
    seq2_lines = FASTAparse(FASTA_file_2)
    scoring_lines = scoring_parse(scoring_file)

    best_score, best_align = local_alignment(seq1_lines, seq2_lines, scoring_lines)

    printNice(best_score, best_align)

if __name__ == "__main__":
    main()
