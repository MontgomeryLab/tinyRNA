""" Script for creating test data for small RNA analysis. """

import numpy as np
import pandas as pd

def parse_org(mature_file, organism="Caenorhabditis elegans", out_file="cel_mature.fa"):
    """ 
    Simply takes the mature miRNA fasta and
    pulls the organism's sequences out
    """
    mature_org = dict()
    with open(mature_file) as mat:
        while True:
            header = mat.readline().strip()
            sequence = mat.readline().strip().replace("U","T")
            
            # EOF
            if not sequence:
                break
            
            if organism in header:
                new_header = header.split()[0]
                mature_org[new_header] = sequence

    with open(out_file, 'w') as out:
        count = 0
        for key, value in mature_org.items():
            # keep just the first 300 for simplicity
            if count < 300:
                out.writelines([key, '\n', value, '\n'])
            else:
                return mature_org

            count +=1

def create_org_table(mature_org,outfile):
    """
    Simply takes the mature miRNA sequences and
    creates a table with the sequences and a count
    to produce.
    """
    np.random.seed(256)
    mature_table = pd.DataFrame.from_dict(list(mature_org.items()))
    mature_table["counts"] = np.random.randint(0,100,size=mature_table.shape[0])
    mature_table.columns = ["name","seq","counts"]

    mature_table.to_csv(outfile)

    return mature_table

def create_fastq(mature_table, outfile, N=100):
    """ 
    Creates a fastq file of N sequences based on
    counts in the sequence table
    """
    i=0
    with open(outfile,'w') as f:
        while i < N:
            for j in range(1,mature_table.loc[i,"counts"]):
                f.writelines('@seq_id_'+str(i)+str(j)+'\n'+mature_table.loc[i,"seq"]+'\n+\n'+'E'*len(mature_table.loc[i,"seq"])+'\n')
            i+=1

def main():
    mature_cel = parse_org(mature_file='testdata/mature.fa')
    mature_table = create_org_table(mature_cel, 'data/mature_counts.csv')
    create_fastq(mature_table, 'cel_mirnas.fq', N=100)

if __name__ == '__main__':
    main()
