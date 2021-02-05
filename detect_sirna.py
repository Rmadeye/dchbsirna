from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt
# import pandas as pd
""" PASENZER/SENSE TO MRNA KURWO JEBANA"""
"""antysense/guide klei sie do riska kurwo jebana"""
"""Colors setup"""
# General colors
CGREEN_BACKGROUND = '\033[102m'
CGREEN = '\033[32m'
CRED = '\033[31m'
CBLUE = '\033[44m'
CEND = '\033[0m'
# Nucleotide colors
nuc_colors = {'A':'\033[32m' + "A" + CEND, 'U':'\033[32m'+'U' +CEND,
              'T':'\033[32m'+'T' +CEND,'G':'\033[31m' + 'G'+'\033[0m','C':'\033[31m'+'C'+'\033[0m' }

class SeqReader:
    def __init__(self, input: str) -> None:
        self.input = input

    def read_fasta_prepare_subsequences(self):
        mRNA_sequence = SeqIO.read(self.input, 'fasta').seq
        subseq_list = [mRNA_sequence[i:i+19] for i in range(0,len(mRNA_sequence)) if len(mRNA_sequence[i:i+19]) == 19]
        return subseq_list

    def create_antisense_strand_from_mRNA_template(self, mrna: Seq) -> Seq:
        return mrna.complement()#+"TT"

    def check_sense_strand_quality(self, mrna: Seq) -> tuple:
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
        siDirect_counter = 0
        reynolds_points = 0
        counter = 0
        requirements_matrix["Sequence"] = mrna._data
        """ AU rich points"""
        if mrna[14:18].count('A') + mrna[14:18].count('U') > 0:
            requirements_matrix["At least 4 A/U on bp 12-19"] = f"{CGREEN}PASSED{CEND}," \
                                                                f" {mrna[14:18].count('A') + mrna[14:18].count('U')}"
            counter += mrna[14:18].count('A') + mrna[14:18].count('U')
            siDirect_counter += mrna[14:18].count('A') + mrna[14:18].count('U')
        """GC content"""
        if (GC(mrna) < 52) and (GC(mrna) > 30):
            siDirect_counter += 1
            reynolds_points += 1
            requirements_matrix["GC content"] = CGREEN + str(int(GC(mrna))) + CEND
        else:
            requirements_matrix["GC content"] = CRED + str(int(GC(mrna))) + CEND

        """Internal repeats"""
        check_reps = 0
        for index, nucleotide in enumerate(mrna):
            try:
                if mrna[index] == mrna[index + 1]:
                    check_reps += 1
            except IndexError:
                break

        if check_reps == 0:
            requirements_matrix["internal repeats"] = CGREEN+'0'+CEND
            siDirect_counter += 1
            reynolds_points += 1
        else:
            requirements_matrix["internal repeats"] = CRED + str(check_reps) + CEND
        """Nucleotides on selected positions"""
        if mrna[2] == 'A':
            counter += 1
            reynolds_points += 1
        if mrna[18] == 'A':
            counter += 1
            reynolds_points += 1
        if mrna[9] == 'U' or mrna[9] == 'T':
            counter += 1
            reynolds_points += 1
        if mrna[12] == 'G':
            counter -= 1

        if (mrna[18] == 'G' or mrna[18] == 'C'):
            counter -= 1

        """ Check GC stretch"""

        gc_stretch_check = [mrna[i:i + 10] for i in range(len(mrna)) if len(mrna[i:i + 10]) == 10]
        for element in gc_stretch_check:
            # print(CRED+element)
            if GC(element) == 0:
                counter -= 1
                requirements_matrix['GC stretch'] = "PRESENT"
                break
            else:
                requirements_matrix['GC stretch'] = CGREEN+"ABSENT"+CEND


        if mrna[0] == "G":
            requirements_matrix["G/C start"] = CGREEN+"PASSED"+CEND
            counter += 1
        else:
            counter = 0

        if mrna[-3] == "A":
            requirements_matrix["A/U on position 19"] = CGREEN+"PASSED"+CEND
            counter += 1
        else:
            counter = 0

        """General melting temperature"""
        requirements_matrix["Melting temp"] = round(mt.Tm_NN(seq = mrna, c_seq=mrna.complement()), 1)
        """3' melting temperature"""
        endat3= mrna[8:18]
        passengerat5=endat3.complement()
        requirements_matrix["3' guide strand melting temp"] = round(mt.Tm_NN(seq = endat3, c_seq=passengerat5),1)
        """5' melting temp"""
        endat5= mrna[0:10]
        passengerat3=endat5.complement()
        requirements_matrix["5' guide strand melting temp"] = round(mt.Tm_NN(seq = endat5, c_seq=passengerat3),1)




        if counter > 1:
            requirements_matrix['Total number of points'] = counter
            # print(requirements_matrix)
            return requirements_matrix, True


    def check_uitei_design(self, sense_strand: Seq) -> tuple:
        uitei_rules = {"A/U not in position 19":"",
                       "G/C on position 1":"",
                       "A/U rich within 10-17":"",
                       "GC stretch < 10":"",
                       "Asymmetrical 5'/3' Tm":"",
                       "Unstable 3' low Tm":"",
                       "Total Tm":"",
                       "Result":""}
        points = 0
        """Rule 1"""
        if (sense_strand[18] != "A") or (sense_strand[18] != "U") or (sense_strand[18] != "T"):
            uitei_rules["A/U not in position 19"] = sense_strand[18]
            points += 1
        else:
            uitei_rules["A/U not in position 19"] = sense_strand[18]
        """Rule 2"""
        if (sense_strand[0] == "G") or (sense_strand[0] == "C"):
            uitei_rules["G/C on position 1"] = sense_strand[0]
            points += 1
        else:
            uitei_rules["G/C on position 1"] = sense_strand[0]
        """Rule 3"""
        if sense_strand[9:16].count('A') + sense_strand[9:16].count('U') + sense_strand[9:16].count('T') > 4:
            uitei_rules['A/U rich within 10-17'] = sense_strand[9:16].count('A') +\
                                                   sense_strand[9:16].count('U') + \
                                                   sense_strand[9:16].count('T')
            points += 1
        else:
            uitei_rules['A/U rich within 10-17'] = sense_strand[9:16].count('A') +\
                                                   sense_strand[9:16].count('U') + \
                                                   sense_strand[9:16].count('T')
        """Rule 4"""
        gc_stretch_check = [sense_strand[i:i + 10] for i in range(len(sense_strand)) if len(sense_strand[i:i + 10]) == 10]
        stretch_count = 0
        for element in gc_stretch_check:
            if GC(element) == 0:
                uitei_rules["GC stretch < 10"] = "Present"
                break
            else:
                stretch_count += 1
        if len(gc_stretch_check) == stretch_count:
            uitei_rules["GC stretch < 10"] = "Absent"
            points += 1
        """Rule 5"""
        tm_at_5 = round(mt.Tm_NN(seq = sense_strand[0:8], c_seq=sense_strand[0:8].complement()),1)
        tm_at_3 = round(mt.Tm_NN(seq = sense_strand[10:18], c_seq=sense_strand[10:18].complement()),1)
        diff = round(abs(tm_at_5-tm_at_3),1)
        if diff > 5:
            uitei_rules["Asymmetrical 5'/3' Tm"] = f"Difference: {diff}C, OK"
            points += 1
        else:
            uitei_rules["Asymmetrical 5'/3' Tm"] = f"Difference: {diff}C, BAD"
        """Rule 6"""
        if tm_at_5-tm_at_3 >= 10:
            uitei_rules["Unstable 3' low Tm"] = f"Passed, Tm at 3' END = {tm_at_3}C, dT > 10"
            points += 1
        else:
            uitei_rules["Unstable 3' low Tm"] = f"Failed, Tm = {tm_at_3}, dT < 10"
        """Rule 7"""
        total_mt = round(mt.Tm_NN(seq = sense_strand, c_seq=sense_strand.complement()),1)
        if total_mt <= 21.5:
            uitei_rules["Total Tm"]=f"{total_mt}C"
            points += 1
        else:
            uitei_rules["Total Tm"] = f"{total_mt}C, higher than 21.5C, FAIL"
        uitei_rules["Result"] = points
        # print(uitei_rules)
        if points > 5:
            return True, uitei_rules
        else:
            return False, uitei_rules

    def check_reynolds_algorithm(self, mrna: Seq) -> tuple:
        reynolds_rules = {"GC content":"",
                          "A/U within 15-19nt":"",
                          "internal repeats":"",
                          "A at 3nt":"",
                          "A at 19nt":"",
                          "U at 10nt":"",
                          "G/C not at 19nt":"",
                          "Result":""}
        return reynolds_rules, True









if __name__ == '__main__':
    testcase = SeqReader("testAG.fasta")
    subseqlist = testcase.read_fasta_prepare_subsequences()
    for mrna in subseqlist:
        antisense_strand = testcase.create_antisense_strand_from_mRNA_template(mrna)
        if testcase.check_uitei_design(mrna)[0]:
            print(f"{CGREEN_BACKGROUND}Sense strand   :  3' {CEND}  {''.join([nuc_colors[x] for x in mrna])}TT  5'")
            print("                       "+"|"*len(antisense_strand))
            print(f"{CBLUE}Antisense strand: 5' {CEND}TT{''.join([nuc_colors[x] for x in antisense_strand])}    3'")
            # for key, value in testcase.check_sense_strand_quality(mrna)[0].items():
            #     print(f'{key}: {value}')
            # print("\n")
            for key, value in testcase.check_uitei_design(mrna)[1].items():
                print(f'{key}: {value}')
            print("\n")


