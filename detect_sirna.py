from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
# import pandas as pd

class SeqReader:
    def __init__(self, input: str) -> None:
        self.input = input

    def read_fasta_prepare_subsequences(self):
        sequence = SeqIO.read(self.input, 'fasta').seq
        subseq_list = [sequence[i:i+21] for i in range(0,len(sequence)) if len(sequence[i:i+21]) == 21]
        return subseq_list

    def create_complimentary_strand(self, mrna: Seq) -> Seq:
        return mrna.complement()+"TT"

    def check_sirna_quality(self, sirna: Seq) -> bool:
        requirements_matrix = {"Sequence":"",
                               "Total number of points":"0",
            "G/C start": "FAILED",
            "A/U on position 19": "FAILED",
            "At least 4 A/U on bp 12-19": "FAILED",
                               "GC content":"0",
                               "internal repeats":"x"
        }
        counter = 0
        requirements_matrix["Sequence"] = sirna._data

        if sirna[14:18].count('A') + sirna[14:18].count('U') > 0:
            requirements_matrix["At least 4 A/U on bp 12-19"] = f"PASSED," \
                                                                f" {sirna[14:18].count('A') + sirna[14:18].count('U')}"
            counter += sirna[14:18].count('A') + sirna[14:18].count('U')

        if (GC(sirna) < 52) and (GC(sirna) > 36):
            counter += 1
            requirements_matrix["GC content"] = int(GC(sirna))

        check_reps = 0
        for index, nucleotide in enumerate(sirna):
            try:
                if sirna[index] == sirna[index+1]:
                    check_reps += 1
            except IndexError:
                break

        if check_reps != 0:
            requirements_matrix["internal repeats"] = check_reps
        else:
            counter += 1

        if sirna[2] == 'A':
            counter += 1
        if sirna[18] == 'A':
            counter += 1
        if sirna[9] == 'U':
            counter += 1
        if sirna[12] == 'G':
            counter -= 1

        if (sirna[18] == 'G' or sirna[18] == 'C'):
            counter -= 1

        if sirna[0] == "G":
            requirements_matrix["G/C start"] = "PASSED"
            counter += 1
        else:
            counter = 0

        if sirna[-3] == "A":
            requirements_matrix["A/U on position 19"] = "PASSED"
            counter += 1
        else:
            counter = 0


        if counter > 7:
            requirements_matrix['Total number of points'] = counter
            print(requirements_matrix)





if __name__ == '__main__':
    testcase = SeqReader("example.fasta")
    subseqlist = testcase.read_fasta_prepare_subsequences()
    for el in subseqlist:
        passenger_strand = testcase.create_complimentary_strand(el)
        if testcase.check_sirna_quality(passenger_strand):
            print(f"Passenger   : 5'{testcase.create_complimentary_strand(el)}  3'")
            print("                |||||||||||||||||||||")
            print(f"Guide strand: 3'{el}    5'"  )
            print("\n")


