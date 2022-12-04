import argparse
import os








if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    #add argument -i:
    parser.add_argument(
        '-i',
        type=str,
        help='input sequence'
    )
    parser.add_argument(
        '-db',
        type=str,
        help='database file (fasta format)',
        default = "tRNAs.fasta"
    )
    parser.add_argument(
        "-E",
        type = int,
        default= 4
        )
    parser.add_argument(
        "-ss",
        type = int
        )
    parser.add_argument(
        "-seed",
        type = str
        )

    args = parser.parse_args()
    query = args.i
    db_path = args.db
    seuil = args.E
    seuil_ss = args.ss
    seed = args.seed

    #======================================================
    # Auxilary functions
    #======================================================

    def get_k_mer(k, seq):
        k_mers = []
        for i in range(len(seq) - k + 1):
            k_mers.append(seq[i:i+k])
        return k_mers

    def read_fasta(fasta_file):
        seqs = [] 
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    continue
                seqs.append(line.strip())
        return seqs

    def exact_research(seed, sequence):
        # find all the positions of the seed in the sequence
        positions = []
        for i in range(len(sequence) - len(seed) + 1):
            if sequence[i:i+len(seed)] == seed:
                positions.append(i)
        return positions

    class Alignment():
        def __init__(
            self, 
            query_align, seq_align, bit_score, e_value,
            query_start, query_end, sequence_start, sequence_end):
            self.query_align = query_align
            self.seq_align = seq_align
            self.bit_score = bit_score
            self.e_value = e_value
            self.query_start = query_start
            self.query_end = query_end
            self.sequence_start = sequence_start
            self.sequence_end = sequence_end
        
        #def get_score(self):
        def __str__(self):
            return f"[query_align: {self.query_align}, \nseq_align: {self.seq_align}, \nbit_score: {self.bit_score}, \ne_value: {self.e_value}, \nquery_start: {self.query_start}, \nquery_end: {self.query_end}, \nsequence_start: {self.sequence_start}, \nsequence_end: {self.sequence_end}]" 




    def extension(alignment, query, sequence):
        
        # extend the alignment to the left and to the right
        # phase 1: extend to the left
        # phase 2: extend to the right

        phases = ["left", "right"]
        for phase in phases:
            current_query_pos = \
            alignment.query_start if phase == "left" else alignment.query_end
            current_sequence_pos = \
            alignment.sequence_start if phase == "left" else alignment.sequence_end
            added_scores = 0
            max_consecutive_mismatches = seuil//4 
            consecutive_mismatches_count = 0
            while True:
                print(current_query_pos)
                direction = -1 if phase == "left" else 1
                if phase == "left":
                    if alignment.query_start == 0 or alignment.sequence_start == 0:
                        break
                if phase == "right":
                    if alignment.query_end == len(query) - 1 or alignment.sequence_end == len(sequence) -1:
                        break

                if query[current_query_pos + direction] == sequence[current_sequence_pos + direction]:
                    added_scores += 5
                    consecutive_mismatches_count = 0
                    current_query_pos += direction
                    current_sequence_pos += direction
                else:
                    consecutive_mismatches_count += 1
                    if consecutive_mismatches_count > max_consecutive_mismatches:
                        break
                    else:
                        added_scores -= 4
                        current_sequence_pos += direction
                        current_query_pos += direction
                
            if added_scores > 0:
                if phase == "left":
                    alignment.query_start = current_query_pos
                    alignment.sequence_start = current_sequence_pos
                if phase == "right":
                    alignment.query_end = current_query_pos
                    alignment.sequence_end = current_sequence_pos

            alignment.query_align = query[alignment.query_start:alignment.query_end + 1]
            alignment.seq_align = sequence[alignment.sequence_start:alignment.sequence_end + 1]
        return alignment



    
    #======================================================
    #======================================================
    #======================================================
    length_seed = len(seed)
    # read the database
    db = read_fasta(db_path)

    def get_alignment(query):

        # get k_mers of the sequence
        k_mers = get_k_mer(length_seed, query)

        # exact research for each sequence of the database
        # HSPs
        HSPs = []
        for sequence in db:
            result = []
            for (index,seed) in enumerate(k_mers):
                positions = exact_research(seed, sequence)
                for position in positions:
                    alignment = Alignment(
                        query_align =  seed,
                        seq_align =  seed,
                        bit_score = 0,
                        e_value = 0,
                        query_start = index,
                        query_end = index + length_seed - 1,
                        sequence_start = position,
                        sequence_end = position + length_seed - 1
                    )
                    alignment = extension(alignment, query, sequence)                
                    result.append(alignment)
            HSPs.append(result)

        for HSP in HSPs:
            for alignment in HSP:
                print(alignment)
    
    querys = read_fasta("unknown.fasta")
    print(querys)
    for query in querys:
        get_alignment(query)
        print("=========================================")
       

        




