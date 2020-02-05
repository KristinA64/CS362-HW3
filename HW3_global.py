import sys


def global_alignment(seq_1, seq_2, score_tuple):
    '''
        Computes the optimal global alignment between sequence 1
        and sequence 2, using the given scoring tuple with the
        following structure:

        score_tuple: [match, mismatch, gap-extend, gap-open]
    '''
    G = [] # match/mismatch table
    E = [] # table where V ends with a gap
    F = [] # table where W ends with a gap
    S = [] # best of all other tables
    alignment_table = [] #"diag", "vert", "hori", "zero"
    #score_tuple: [match, mismatch, gap-extend, gap-open]

    seq_1_padded = "-" + seq_1
    seq_2_padded = "-" + seq_2

    #Makes sure that squares with this value are never preferable
    min_score = (score_tuple[3]+score_tuple[2])*len(seq_1)*len(seq_2)

    for x in range(len(seq_1)+1):
        G.append([min_score]*(len(seq_2)+1))
        E.append([min_score]*(len(seq_2)+1))
        F.append([min_score]*(len(seq_2)+1))
        S.append([min_score]*(len(seq_2)+1))
        alignment_table.append(["zero"]*(len(seq_2)+1))

    S[0][0] = 0
    G[0][0] = 0

    for i in range(len(seq_1)+1):
        F[i][0] = score_tuple[3]+(i*score_tuple[2])
        if i>=1:
            S[i][0] = F[i][0]
            alignment_table[i][0] = "vert"

    for j in range(len(seq_2)+1):
        E[0][j] = score_tuple[3] + (j*score_tuple[2])
        if j>=1:
            S[0][j] = E[0][j]
            alignment_table[0][j] = "hori"

    for i in range(1,len(seq_1)+1):
        for j in range(1,len(seq_2)+1):
            recur_helper(G, E, F, S, i, j, seq_1_padded[i], seq_2_padded[j], score_tuple, alignment_table)

    best_pos =  (len(seq_1),len(seq_2))
    best_score = S[len(seq_1)][len(seq_2)]

    best_align = TraceBack(best_pos, alignment_table, seq_1_padded, seq_2_padded)



    return best_score, best_align

def recur_helper(G, E, F, S, i, j, cur_char_1, cur_char_2, score_tuple, alignment_table):
    '''
        Performs the table filling for a single cell in all 5 tables.
    '''
    match_or_mismatch = 0 #assume match
    if cur_char_1 != cur_char_2:
        match_or_mismatch = 1 #it's a mismatch

    ### G ###
    G[i][j] = S[i-1][j-1] + score_tuple[match_or_mismatch]

    ### E ###
    E[i][j] = max(E[i][j-1]+score_tuple[2], S[i][j-1]+score_tuple[3]+score_tuple[2])

    ### F ###
    F[i][j] = max(F[i-1][j]+score_tuple[2], S[i-1][j]+score_tuple[3]+score_tuple[2])

    ### S + Alignment Table ###
    diag_score = G[i][j]
    vert_score = F[i][j]
    hori_score = E[i][j]

    if hori_score >= diag_score and hori_score >= vert_score:
        alignment_table[i][j] = "hori"
        S[i][j] = hori_score
    elif diag_score >= vert_score and diag_score >= hori_score:
        alignment_table[i][j] = "diag"
        S[i][j] = diag_score
    elif vert_score >= diag_score and vert_score >= hori_score:
        alignment_table[i][j] = "vert"
        S[i][j] = vert_score

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
        return [int(i) for i in score_tuple] #score_tuple: [match, mismatch, gap-extend, gap-open]

def printNice(best_score, best_align):
    '''
        Takes an optimal score and a best alignment pair,
        and prints them in a readable fashion.
    '''
    print("Alignment score:", best_score)
    print("Optimal Global Alignment with Affine Gap:")
    print(best_align[0])
    print(best_align[1])

def main():
    FASTA_file_1 = sys.argv[1]
    FASTA_file_2 = sys.argv[2]
    scoring_file = sys.argv[3]

    seq1_lines = FASTAparse(FASTA_file_1)
    seq2_lines = FASTAparse(FASTA_file_2)
    scoring_lines = scoring_parse(scoring_file)

    best_score, best_align = global_alignment(seq1_lines, seq2_lines, scoring_lines)

    printNice(best_score, best_align)
if __name__ == "__main__":
    main()
