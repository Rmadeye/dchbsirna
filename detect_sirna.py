from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt
import termcolor
# import pandas as pd

"""Colors setup"""
# General colors
CGREEN_BACKGROUND = '\033[102m'
CGREEN = '\033[32m'
CRED = '\033[31m'
CBLUE = '\033[44m'
CEND = '\033[0m'
# Nucleotide colors
nuc_colors = {'A':'\033[32m' + "A" + CEND,
              'T':'\033[32m'+'T' +CEND,'G':'\033[31m' + 'G'+'\033[0m','C':'\033[31m'+'C'+'\033[0m' }

class SeqReader:
    def __init__(self, input: str) -> None:
        self.input = input

    def read_fasta_prepare_subsequences(self):
        sequence = SeqIO.read(self.input, 'fasta').seq
        subseq_list = [sequence[i:i+21] for i in range(0,len(sequence)) if len(sequence[i:i+21]) == 21]
        return subseq_list

    def create_complimentary_strand(self, mrna: Seq) -> Seq:
        return mrna.complement()#+"TT"

    def check_sirna_quality(self, sirna: Seq) -> tuple:
        requirements_matrix = {"Sequence":"",
                               "Total number of points":"",
            "G/C start": "",
            "A/U on position 19": "FAILED",
            "At least 4 A/U on bp 12-19": "FAILED",
                               "GC content":"0",
                               "internal repeats":"",
                               "GC stretch": "absent",
                               "Melting temp":"",
                               "3' guide strand melting temp":"",
                               "5' guide strand melting temp": ""

        }
        counter = 0
        requirements_matrix["Sequence"] = sirna._data

        if sirna[14:18].count('A') + sirna[14:18].count('U') > 0:
            requirements_matrix["At least 4 A/U on bp 12-19"] = f"{CGREEN}PASSED{CEND}," \
                                                                f" {sirna[14:18].count('A') + sirna[14:18].count('U')}"
            counter += sirna[14:18].count('A') + sirna[14:18].count('U')

        if (GC(sirna) < 52) and (GC(sirna) > 30):
            counter += 1
            requirements_matrix["GC content"] = CGREEN+str(int(GC(sirna)))+CEND
        else:
            requirements_matrix["GC content"] = CRED+str(int(GC(sirna)))+CEND


        check_reps = 0
        for index, nucleotide in enumerate(sirna):
            try:
                if sirna[index] == sirna[index+1]:
                    check_reps += 1
            except IndexError:
                break

        if check_reps == 0:
            requirements_matrix["internal repeats"] = CGREEN+'0'+CEND
            counter += 1
        else:
            requirements_matrix["internal repeats"] = CRED + str(check_reps) + CEND
            counter -= 1

        if sirna[2] == 'A':
            counter += 1
        if sirna[9] == 'U' or sirna[9] == 'T':
            counter += 1
        if sirna[18] == 'A':
            counter += 1
        if sirna[9] == 'U' or sirna[9] == 'T':
            counter += 1
        if sirna[12] == 'G':
            counter -= 1

        if (sirna[18] == 'G' or sirna[18] == 'C'):
            counter -= 1

        """ Check GC stretch"""

        gc_stretch_check = [sirna[i:i+9] for i in range(len(sirna)) if len(sirna[i:i+9]) == 9]
        for element in gc_stretch_check:
            if GC(element) == 0:
                counter -= 1
                requirements_matrix['GC stretch'] = "PRESENT"
                break
            else:
                requirements_matrix['GC stretch'] = CGREEN+"ABSENT"+CEND


        if sirna[0] == "G":
            requirements_matrix["G/C start"] = CGREEN+"PASSED"+CEND
            counter += 1
        else:
            counter = 0

        if sirna[-3] == "A":
            requirements_matrix["A/U on position 19"] = CGREEN+"PASSED"+CEND
            counter += 1
        else:
            counter = 0

        """General melting temperature"""
        requirements_matrix["Melting temp"] = round(mt.Tm_NN(seq = sirna, c_seq=sirna.complement()),1)
        """3' melting temperature"""
        endat3=sirna[8:18]
        passengerat5=endat3.complement()
        requirements_matrix["3' guide strand melting temp"] = round(mt.Tm_NN(seq = endat3, c_seq=passengerat5),1)
        # print(endat3+'TT')
        # print("||||||")
        # print(passengerat5)
        """5' melting temp"""
        endat5=sirna[0:10]
        passengerat3=endat5.complement()
        requirements_matrix["5' guide strand melting temp"] = round(mt.Tm_NN(seq = endat5, c_seq=passengerat3),1)




        if counter > 6:
            requirements_matrix['Total number of points'] = counter
            # print(requirements_matrix)
            return requirements_matrix, True





if __name__ == '__main__':
    testcase = SeqReader("example.fasta")
    subseqlist = testcase.read_fasta_prepare_subsequences()
    for el in subseqlist:
        passenger_strand = testcase.create_complimentary_strand(el)
        if testcase.check_sirna_quality(el):
            print(f"{CGREEN_BACKGROUND}Passenger   : 3' {CEND}  {''.join([nuc_colors[x] for x in testcase.create_complimentary_strand(el)])}TT  5'")
            print("                   |||||||||||||||||||||")
            print(f"{CBLUE}Guide strand: 5' {CEND}TT{''.join([nuc_colors[x] for x in el])}    3'")
            for key, value in testcase.check_sirna_quality(el)[0].items():
                print(f'{key}: {value}')
            print("\n")


