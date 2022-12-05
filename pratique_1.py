import os

def read_fasta(fasta_file):
    seqs = [] 
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            seqs.append(line.strip())
    return seqs

if __name__ == '__main__':
    querys = read_fasta("unknown.fasta")
    for query in querys:
       os.system(f"python3 plast.py -i {query} -db tRNAs.fasta -E 4 -ss 0.1 -seed '11111111111'")
