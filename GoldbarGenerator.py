""""
Generating Goldbar for Genetic Circuits with certain Rules

(P) Promoter Road Blocking: Road-blocking promoters must be alone or infront of regular promoters
(L) Leaky Terminators: Leaky Terminators can only be at end of circuit
(O) Orthogonality: Certain parts can not be in the same circuit
(I) Part-junction Interference: Certain parts can not be next to each other
(R) Do not Reuse Parts: Parts will not be reused in the same sequence

Inputs:
Principles to include
Number of sequences
Number of Cycles
Part Library - txt file
Not Orthogonal Library - csv file

Outputs:
Txt file of GOLDBAR
Txt file of Categories
"""

import numpy as np
import pandas as pd
import csv
import json

FILEPATH = input('Enter path to Part Library: ')
part_library = pd.read_csv(FILEPATH, sep=',')
print(part_library)

number_of_tus = int(input('\nEnter number of transcriptional units for each sequence: '))

principles = input('\nEnter principles to use separated by a space [P L O I R]: '
                   '\n(P) Promoter Road Blocking: Road-blocking promoters must be alone or infront of regular promoters '
                   '\n(L) Leaky Terminators: Leaky Terminators can only be at end of circuit '
                   '\n(O) Not Orthogonal: Certain parts can not be in the same circuit '
                   '\n(I) Part-junction Interference: Certain parts can not be next to each other '
                   '\n(R) Do not Reuse Parts: Parts will not be reused in the same sequence'
                   '\nEnter: ')


def principles_used(principles):
    """
    Takes the principles string and returns True or False for each principle. Also returns orthogonal library if
    Not Orthogonal (O) is true.

    :param principles: string
    :return : P, L, O, I, R, orthogonal_library, pji_library, repeated_library
    """

    # Turn on principles
    principles = principles.lower()
    principles = principles.split()

    if 'p' in principles:
        P = True
    else:
        P = False
    if 'l' in principles:
        L = True
    else:
        L = False
    if 'o' in principles:
        O = True
        FILEPATH = input('\nEnter path to Not Orthogonal Library: ')
        orthogonal_library = pd.read_csv(FILEPATH, sep=',')
        print("\nNot Orthogonal Library:")
        print(orthogonal_library)
    else:
        O = False
        orthogonal_library = None
    if 'i' in principles:
        I = True
        FILEPATH = input('\nEnter path to Part Junction Interference Library: ')
        pji_library = pd.read_csv(FILEPATH, sep=',')
        print("\nPart Junction Interference Library:")
        print(pji_library)
    else:
        I = False
        pji_library = None
    if 'r' in principles:
        R = True
    return P, L, O, I, R, orthogonal_library, pji_library


def extract_parts(part_library):
    """
    Extracts parts from part_library dataframe and returns the library as a dictionary.

    :param part_library:
    :return : parts
    """
    # Promoters
    blocking_promoters = part_library['promoter RB'].dropna().to_list()
    regular_promoters = part_library['promoter'].dropna().to_list()
    all_promoters = np.concatenate((regular_promoters, blocking_promoters)).tolist()

    # RBS
    rbs = part_library['rbs'].dropna().to_list()

    # CDS
    cds = part_library['cds'].dropna().to_list()

    # Terminators
    leaky_terminators = part_library['L terminators'].dropna().to_list()
    regular_terminators = part_library['terminators'].dropna().to_list()
    all_terminators = np.concatenate((regular_terminators, leaky_terminators)).tolist()

    # Parts Not to be repeated
    not_repeated_parts = part_library['DNR'].dropna().to_list()

    # All Parts
    all_parts = all_promoters + rbs + cds + all_terminators

    parts = {'all': all_parts,
             'promoter': all_promoters,
             'ribosomeBindingSite': rbs,
             'cds': cds,
             'terminator': all_terminators,
             'promoter_notroadblocking': regular_promoters,
             'terminator_notleaky': regular_terminators,
             'do not repeat': not_repeated_parts}

    return parts


def goldbar_generator(principles, number_of_tus, part_library):
    """
    Generates Goldbar and Categories. Creates txt file for Goldbar and
    Categories.

    :param principles:
    :param number_of_tus:
    :param part_library:
    :return: Nothing
    """

    # Extract Parts
    parts = extract_parts(part_library)

    # activate principles
    P, L, O, I, R, orthogonal_library, pji_library = principles_used(principles)

    goldbar = []
    categories = {"promoter_any": {"promoter": parts["promoter"]},
                  "rbs": {"ribosomeBindingSite": parts["ribosomeBindingSite"]},
                  "cds": {"cds": parts["cds"]},
                  "terminator_any": {"terminator": parts["terminator"]}}

    # Concrete number of transcription units
    seq = []
    for i in range(number_of_tus):

        seq.append("promoter_any then zero-or-one(promoter_any) then one-or-more(rbs then cds) then terminator_any")

    seq = ' then '.join(seq)
    seq = f"({seq})"
    goldbar.append(seq)


    if L:
        # Add to goldbar
        goldbar.append("(zero-or-more(any_except_terminator_leaky) then zero-or-one(terminator_any))")

        # add any_except_terminator_notleaky to categories
        categories["any_except_terminator_leaky"] = {"promoter": parts["promoter"],
                                                        "ribosomeBindingSite": parts["ribosomeBindingSite"],
                                                        "cds": parts["cds"],
                                                        "terminator": parts["terminator_notleaky"]}

    if P:
        # Add to goldbar
        goldbar.append("(one-or-more(promoter_any then zero-or-one(promoter_notroadblocking) then rbs then one-or-more(cds) then terminator_any))")

        # Add promoter_notroadblocking to categories
        categories["promoter_notroadblocking"] = {"promoter": parts["promoter_notroadblocking"]}

    if O:
        # Extract orthogonal dictionary
        orthogonal_library_dict = orthogonal_library.to_dict('list')


        for part in list(orthogonal_library_dict.keys()):

            # Add part to categories
            if part in parts["promoter"]:
                categories[part] = {"promoter": [part]}

            elif part in parts["ribosomeBindingSite"]:
                categories[part] = {"ribosomeBindingSite": [part]}

            elif part in parts["cds"]:
                categories[part] = {"cds": [part]}

            else:
                categories[part] = {"terminator": [part]}

            # Add any_except_part to categories
            categories[f"any_except_{part}"] = {"promoter": [x for x in parts["promoter"] if x != part],
                                                 "ribosomeBindingSite": [x for x in parts["ribosomeBindingSite"] if x != part],
                                                 "cds": [x for x in parts["cds"] if x != part],
                                                 "terminator": [x for x in parts["terminator"] if x != part]}

            for part1 in orthogonal_library[part].dropna().to_list():
                # Add to goldbar
                goldbar.append(f"zero-or-more(any_except_{part}and{part1}) then zero-or-one(({part} then zero-or-more(any_except_{part1})) or ({part1} then zero-or-more(any_except_{part})))")

                # Add part1 to categories
                if part1 in parts["promoter"]:
                    categories[part1] = {"promoter": [part1]}
                    except_part1 = parts

                elif part1 in parts["ribosomeBindingSite"]:
                    categories[part1] = {"ribosomeBindingSite": [part1]}

                elif part1 in parts["cds"]:
                    categories[part1] = {"cds": [part1]}

                else:
                    categories[part1] = {"terminator": [part1]}

                # Add any_except_part1 to categories
                categories[f"any_except_{part1}"] = {"promoter": [x for x in parts["promoter"] if x != part1],
                                                     "ribosomeBindingSite": [x for x in parts["ribosomeBindingSite"] if x != part1],
                                                     "cds": [x for x in parts["cds"] if x != part1],
                                                     "terminator": [x for x in parts["terminator"] if x != part1]}

                # Add any_except_part+part1 to categories
                categories[f"any_except_{part}and{part1}"] = {"promoter": [x for x in parts["promoter"] if x not in [part, part1]],
                                                     "ribosomeBindingSite": [x for x in parts["ribosomeBindingSite"] if x not in [part, part1]],
                                                     "cds": [x for x in parts["cds"] if x not in [part, part1]],
                                                     "terminator": [x for x in parts["terminator"] if x not in [part, part1]]}
    if I:
        # Extract pji dictionary
        pji_library_dict = pji_library.to_dict('list')


        for part in list(pji_library_dict.keys()):

            # Add part to categories
            if part in parts["promoter"]:
                categories[part] = {"promoter": [part]}

            elif part in parts["ribosomeBindingSite"]:
                categories[part] = {"ribosomeBindingSite": [part]}

            elif part in parts["cds"]:
                categories[part] = {"cds": [part]}

            else:
                categories[part] = {"terminator": [part]}

            # Add any_except_part to categories
            categories[f"any_except_{part}"] = {"promoter": [x for x in parts["promoter"] if x != part],
                                                 "ribosomeBindingSite": [x for x in parts["ribosomeBindingSite"] if x != part],
                                                 "cds": [x for x in parts["cds"] if x != part],
                                                 "terminator": [x for x in parts["terminator"] if x != part]}

            for part1 in pji_library[part].dropna().to_list():
                # Add to goldbar
                goldbar.append(f"((one-or-more(any_except_{part}and{part1})) or (one-or-more(zero-or-more(any_except_{part}and{part1}) then any_except_{part1} then {part} then any_except_{part1} then zero-or-more(any_except_{part}and{part1}))) or (one-or-more(zero-or-more(any_except_{part}and{part1}) then any_except_{part} then {part1} then any_except_{part} then zero-or-more(any_except_{part}and{part1}))))")

                # Add part1 to categories
                if part1 in parts["promoter"]:
                    categories[part1] = {"promoter": [part1]}
                    except_part1 = parts

                elif part1 in parts["ribosomeBindingSite"]:
                    categories[part1] = {"ribosomeBindingSite": [part1]}

                elif part1 in parts["cds"]:
                    categories[part1] = {"cds": [part1]}

                else:
                    categories[part1] = {"terminator": [part1]}

                # Add any_except_part1 to categories
                categories[f"any_except_{part1}"] = {"promoter": [x for x in parts["promoter"] if x != part1],
                                                     "ribosomeBindingSite": [x for x in parts["ribosomeBindingSite"] if x != part1],
                                                     "cds": [x for x in parts["cds"] if x != part1],
                                                     "terminator": [x for x in parts["terminator"] if x != part1]}

                # Add any_except_part+part1 to categories
                categories[f"any_except_{part}and{part1}"] = {"promoter": [x for x in parts["promoter"] if x not in [part, part1]],
                                                     "ribosomeBindingSite": [x for x in parts["ribosomeBindingSite"] if x not in [part, part1]],
                                                     "cds": [x for x in parts["cds"] if x not in [part, part1]],
                                                     "terminator": [x for x in parts["terminator"] if x not in [part, part1]]}
                
    if R:
        for part in parts['do not repeat']:
            # Add to goldbar
            goldbar.append(f"(zero-or-more(any_except_{part}) then zero-or-one({part} then zero-or-more(any_except_{part})))")

            # Add part to categories
            if part in parts["promoter"]:
                categories[part] = {"promoter": [part]}

            elif part in parts["ribosomeBindingSite"]:
                categories[part] = {"ribosomeBindingSite": [part]}

            elif part in parts["cds"]:
                categories[part] = {"cds": [part]}

            else:
                categories[part] = {"terminator": [part]}

            categories[f"any_except_{part}"] = {"promoter": [x for x in parts["promoter"] if x != part],
                                                 "ribosomeBindingSite": [x for x in parts["ribosomeBindingSite"] if x != part],
                                                 "cds": [x for x in parts["cds"] if x != part],
                                                 "terminator": [x for x in parts["terminator"] if x != part]}


    # Join goldbar together with "and0"
    goldbar = " and0 ".join(goldbar)

    # print goldbar and categories
    print(f'\nGoldbar:\n{goldbar}')
    print(f'\nCategories:')
    categories_json = json.dumps(categories)
    print(categories_json)

    write_txt = input('\nOutput Goldbar and Categories as txt file? (y/n): ')

    if write_txt.lower() == 'y':
        filename = input('Enter name for file (without .txt): ') + '.txt'
        print(f'Filename: {filename}')
        f = open(filename, "w")

        f.write("Goldbar:\n")
        f.write(goldbar)
        f.write('\n\n')

        goldbar_split = goldbar.split("and1 ")
        f.write("Goldbar (rule by rule):\n")
        for i in range(len(goldbar_split)):
            f.write(goldbar_split[i])
            if i != len(goldbar_split) - 1:
                f.write('and1\n')

        f.write("\n\n")
        categories_split = categories_json.split('},')
        f.write("Categories:\n")
        for i in range(len(categories_split)):
            f.write(categories_split[i])
            if i != len(categories_split) - 1:
                f.write('},\n')
        f.close()
        print('\nFile Created!')


goldbar_generator(principles, number_of_tus, part_library)
