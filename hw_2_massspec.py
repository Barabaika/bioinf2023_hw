import sys
import argparse

'''
Script for predicting cyclic peptide from txt file with tab separated experiment spectrum values.
No special env needed.

ARGS:
1- input file

KWARGS:
-n --top_n (default = 10) : number of top peptides in leaderboard to save at each iteration
-c --cyclic               : predict cyclic peptide or not

NOTES:
- spectrum values are rounded to integers
- no output file, only print in consol  

USAGE:
python hw_2_massspec.py data/Spectrum_task2.txt -n 10 --cyclic
'''

letters_mass = {
    'A':	71.03711,
    'R':	156.10111,
    'N':	114.04293,
    'D':	115.02694,
    'C':	103.00919,
    'E':	129.04259,
    'Q':	128.05858,
    'G':	57.02146,
    'H':	137.05891,
    'I':	113.08406,
    'L':	113.08406,
    'K':	128.09496,
    'M':	131.04049,
    'F':	147.06841,
    'P':	97.05276,
    'S':	87.03203,
    'T':	101.04768,
    'W':	186.07931,
    'Y':	163.06333,
    'V':	99.06841}


def mass(
        peptide: str or list
        ):
    
    return round(sum([letters_mass[l] for l in peptide]), 0)

def expand(liderboard):
    expanded_liderboard = []
    for peptide in liderboard:
        for letter in letters_mass:
            expanded_liderboard.append(peptide + letter)
    new_max_length = len(peptide + letter)

    return expanded_liderboard, new_max_length


def score_cyclic(
        peptide: str or list, 
        spectrum: dict
        ):

    score = 0
    tmp_peptide = peptide*2
    tmp_spectrum = spectrum.copy()
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)+ i+1):
            substr = tmp_peptide[i:j]
            sub_mass = mass(substr)
            if sub_mass in tmp_spectrum:
                score +=1
                tmp_spectrum[sub_mass] -=1
                if tmp_spectrum[sub_mass] <= 0:
                    del tmp_spectrum[sub_mass]

    return score


def score_uncyclic(
        peptide: str or list, 
        spectrum: dict
        ):

    score = 0
    tmp_peptide = peptide
    tmp_spectrum = spectrum.copy()
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)):
            substr = tmp_peptide[i:j]
            sub_mass = mass(substr)
            if sub_mass in tmp_spectrum:
                score +=1
                tmp_spectrum[sub_mass] -=1
                if tmp_spectrum[sub_mass] <= 0:
                    del tmp_spectrum[sub_mass]

    return score


def trim(
        liderboard: list, 
        spectrum: dict, 
        N: int
        ):
    scores = []
    for_print = []
    for j in range(1, len(liderboard)):
        peptide = liderboard[j]
        scores.append((score(peptide, spectrum), j))
        for_print.append((peptide, score(peptide, spectrum)))

    scores_sorted = sorted(scores, reverse=True)
    liderboard_sorted = [liderboard[score_tuple[1]] for score_tuple in scores_sorted]
    
    for j in range(N+1, len(liderboard_sorted)):
        
        if scores_sorted[j][0] < scores_sorted[N][0]:
            liderboard_sorted = liderboard_sorted[: j]
            return liderboard_sorted
    return liderboard_sorted

def liderboard_peptide_seqiencing(
        spectrum: dict, 
        N: int,
        score
        ):
    liderboard = ['']
    leader_peptide = ''
    leader_score = score(leader_peptide, spectrum)
    continue_flag = True
    while len(liderboard) != 0:
        liderboard, _ = expand(liderboard)
        tmp_leaderboard = liderboard.copy()
        for peptide in tmp_leaderboard:
            m_peptide = mass(peptide)
            if m_peptide == max(spectrum):
                tmp_s = score(peptide, spectrum)
                if tmp_s > leader_score:
                    leader_peptide = peptide
                    leader_score = tmp_s
                    print('new leader peptide', leader_peptide)
                    print(f'new leader score: {leader_score}/{sum(spectrum.values())}')
            elif m_peptide > max(spectrum):
                liderboard.remove(peptide)

        liderboard = trim(liderboard, spectrum, N)
    print('\nfinal result:')
    print('PEPTIDE', leader_peptide)
    print(f'SCORE: {leader_score}/{sum(spectrum.values())}')
    return leader_peptide, leader_score
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument('-n', '--top_n', default=10)
    parser.add_argument('-c', '--cyclic', action="store_true")
    args = parser.parse_args()

    print(f'start process {args.input_file}\n')

    # define score function
    score = score_uncyclic
    if args.cyclic:
        score = score_cyclic

    # load input file to list
    spectr_list = []
    with open(args.input_file, 'r') as f:
        for line in f:
            line_split = line.split(' ')
            line_split = [float(i) for i in line_split]
            spectr_list += line_split
    # process input list to input dict
    spectr_list = [round(peak, 0) for peak in spectr_list]
    spectr_dict = {peak:0 for peak in spectr_list}
    for k in spectr_list:
        spectr_dict[k] += 1

    # run algoritm
    res_pep, res_score = liderboard_peptide_seqiencing(spectr_dict, int(args.top_n), score)


    

    
    
