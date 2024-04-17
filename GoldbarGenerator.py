""""
Generating Goldbar for Genetic Circuits with certain Rules

(P) Promoter Road Blocking: Road-blocking promoters must be alone or infront of regular promoters
(L) Leaky Terminators: Leaky Terminators can only be at end of circuit
(O) Orthogonality: Certain parts can not be in the same circuit
(I) Part-junction Interference: Certain parts can not be next to each other 

(R) Do not Reuse Certain Parts: Certain Parts will not be reused in the same sequence

(T) Certain Parts Must be in the same circuit (Needs Work)
(M) Parts that have to be included
(E) Place Part at End

Inputs:
Part Library : CSV File - Required
Number of Transcription Units : Int - Optional
Principles : String - Optional (Default is None)
Not Orthogonal Library : CSV File - Optional 
Part Junction Interference Library : CSV File
Together Library : CSV File
Genetic Circuit Format : String - Optional
Any Part Format : String - Optional

Outputs:
Txt file of GOLDBAR
Txt file of Categories
"""

import numpy as np
import pandas as pd
import csv
import json

def principles_used(principles):
    """
    Takes the principles string and returns True or False for each principle. Also returns orthogonal library if
    Not Orthogonal (O) is true.

    :param principles: string
    :return : P, L, O, I, R, T, M, E, orthogonal_library, pji_library, together_library
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
    else:
        R = False

    if 't' in principles:
        T = True
        FILEPATH = input('\nEnter path to Together Library: ')
        together_library = pd.read_csv(FILEPATH, sep=',')
        print("\nTogether Library:")
        print(together_library)
    else:
        T = False
        together_library = None

    if 'm' in principles:
        M = True
    else:
        M = False
    
    if 'e' in principles:
        E = True
    else:
        E = False

    return P, L, O, I, R, T, M, E, orthogonal_library, pji_library, together_library

def extract_parts(part_library):
    """
    Extracts parts from part_library dataframe and returns the library as a dictionary.

    :param part_library: DataFrame
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

    # Parts that have to be in the design
    must_include_parts = part_library['Must_Include'].dropna().to_list()

    # Part that have to be at the end of a design
    end_part = part_library['END'].dropna().to_list()[0]

    # All Parts
    all_parts = all_promoters + rbs + cds + all_terminators

    parts = {'all': all_parts,
             'promoter': all_promoters,
             'ribosomeBindingSite': rbs,
             'cds': cds,
             'terminator': all_terminators,
             'promoter_notroadblocking': regular_promoters,
             'terminator_notleaky': regular_terminators,
             'terminator_leaky': leaky_terminators,
             'do not repeat': not_repeated_parts,
             'must include': must_include_parts,
             'end_part' : end_part}

    return parts

def initialize_categories(parts, allIndividual):
    """
    Initializes Categories Dictionary

    :param parts: Dictionary of Part Types with List of Part Names
    :param allIndividual: String (y/n)
    :return: Dictionary of Categories
    """
    part_types = ["promoter", "ribosomeBindingSite", "cds", "terminator"]

    # Initialize Categories with individual part types
    categories = {"promoter_any": {"promoter": parts["promoter"]},
                  "rbs": {"ribosomeBindingSite": parts["ribosomeBindingSite"]},
                  "cds": {"cds": parts["cds"]},
                  "terminator_any": {"terminator": parts["terminator"]}}
    
    # Add any_part Category (Abstract)
    categories['any_part_abstract'] = {"promoter": [],
                              "ribosomeBindingSite": [],
                              "cds": [],
                              "terminator": []}
    
    # Add any_part_concrete Category (Concrete)
    categories['any_part_concrete'] = {}
    for part_type in part_types:
        if parts[part_type]:
            categories['any_part_concrete'][part_type] = parts[part_type]
    
    if 'n' not in allIndividual.lower():
        for part_type in part_types:
            for part in parts[part_type]:
                categories[part] = {part_type : [part]}
    
    return categories

def add_part(part, parts, categories):
    """
    Adds part to categories

    :param part: String or List, Part to add to categories
    :param parts: Dictionary of Part Types with List of Part Names
    :param categories: Dictionary of Categories
    :return: Dictionary of Updated Categories
    """

    if type(part) == str:
        if part in parts["promoter"]:
            categories[part] = {"promoter": [part]}

        elif part in parts["ribosomeBindingSite"]:
            categories[part] = {"ribosomeBindingSite": [part]}

        elif part in parts["cds"]:
            categories[part] = {"cds": [part]}

        else:
            categories[part] = {"terminator": [part]}
    
    else:
        name = "and".join(part)
        categories[name] = {}

        for p in part:
            if p in parts["promoter"]:
                if categories[name].get("promoter", False):
                    categories[name]["promoter"].append(p)
                else:
                    categories[name]["promoter"] = [p]

            elif p in parts["ribosomeBindingSite"]:
                if categories[name].get("ribosomeBindingSite", False):
                    categories[name]["ribosomeBindingSite"].append(p)
                else:
                    categories[name]["ribosomeBindingSite"] = [p]

            elif p in parts["cds"]:
                if categories[name].get("cds", False):
                    categories[name]["cds"].append(p)
                else:
                    categories[name]["cds"] = [p]  

            elif p in parts["terminator"]:
                if categories[name].get("terminator", False):
                    categories[name]["terminator"].append(p)
                else:
                    categories[name]["terminator"] = [p]

    return categories

def any_except_part(part, parts, categories, isConcrete):
    """
    Adds any_except_part to categories

    :param part: Part to add to categories
    :param parts: Dictionary of Part Types with List of Part Names
    :param categories: Dictionary of Categories
    :param isConcrete: Boolean
    :return: Dictionary of Updated Categories
    """

    part_types = ["promoter", "ribosomeBindingSite", "cds", "terminator"]
    
    categories[f"any_except_{part}"] = {}

    for part_type in part_types:
        if (parts[part_type]) or (not isConcrete):
            categories[f"any_except_{part}"][part_type] = [x for x in parts[part_type] if x != part]

    return categories

def any_except_part1andpart2(part1, part2, parts, categories, isConcrete):
    """
    Adds any_except_part1andpart2 to categories

    :param part: 1st Part to add to categories
    :param part: 2nd Part to add to categories
    :param parts: Dictionary of Part Types with List of Part Names
    :param categories: Dictionary of Categories
    :param isConcrete: Boolean
    :return: Dictionary of Updated Categories
    """

    part_types = ["promoter", "ribosomeBindingSite", "cds", "terminator"]

    categories[f"any_except_{part1}and{part2}"] = {}

    for part_type in part_types:
        if (parts[part_type]) or (not isConcrete):
            categories[f"any_except_{part1}and{part2}"][part_type] = [x for x in parts[part_type] if x not in [part1, part2]]

    return categories

def any_except_multiple_parts(except_parts, parts, categories, isConcrete):
    """
    Adds any_except_parts to categories

    :param except_parts: List of parts to add to categories
    :param parts: Dictionary of Part Types with List of Part Names
    :param categories: Dictionary of Categories
    :param isConcrete: Boolean
    :return: Dictionary of Updated Categories
    """

    part_types = ["promoter", "ribosomeBindingSite", "cds", "terminator"]

    name = "any_except_" + "and".join(except_parts)
    categories[name] = {}

    for part_type in part_types:
        if (parts[part_type]) or (not isConcrete):
            categories[name][part_type] = [x for x in parts[part_type] if x not in except_parts]

    return categories


def goldbar_generator(principles, number_of_tus, part_library, circuit_format, any_part_format, allIndividual):
    """
    Generates Goldbar and Categories. Creates txt file for Goldbar and
    Categories.

    :param principles:
    :param number_of_tus:
    :param part_library:
    :param circuit_format:
    :param any_part_format: string
    :return: Nothing
    """

    # Initialize GOLDBAR Dictionary
    goldbar = {'N' : [],
               'L' : [],
               'P' : [],
               'O' : [],
               'I' : [],
               'R' : [],
               'T' : [],
               'M' : [],
               'E' : []}

    # Extract Parts
    parts = extract_parts(part_library)

    # activate principles
    P, L, O, I, R, T, M, E, orthogonal_library, pji_library, together_library = principles_used(principles)

    # Concrete or Abstract any_part
    if 'concrete' in any_part_format.lower():
        any_part_format = 'any_part_concrete'
        isConcrete = True
    else:
        any_part_format = 'any_part_abstract'
        isConcrete = False

    # Initialize Categories
    categories = initialize_categories(parts, allIndividual)

    if number_of_tus:
        # Concrete number of transcription units
        if 'cello' in circuit_format:

            num_of_tus = ''
            greater_than = False
            less_than = False

            for char in number_of_tus:
                if char == '<':
                    less_than = True
                elif char == '>':
                    greater_than = True
                elif char.isdigit():
                    num_of_tus += char

            seq = []
            next = ''
            for i in range(int(num_of_tus)):            
                if less_than:
                    if i == 0:
                        next += "(any_part then any_part) then "
                    else:
                        next += "(any_part then any_part then any_part) then "
                    seq.append(next)
                elif greater_than:
                    seq.append("one-or-more(promoter_any then cds then zero-or-one(terminator_any))")
                else:
                    if i == 0:
                        seq.append("(any_part then any_part)")
                    else:
                        seq.append("(any_part then any_part then any_part)")

            if less_than:
                for i in range(len(seq)):
                    seq[i] = seq[i][:-6]
                    seq[i] = f"({seq[i]})"
                seq = ' or '.join(seq)
            
            else:
                seq = ' then '.join(seq)

            seq = f"({seq})"
            goldbar['N'].append(seq)

        elif 'pc' in circuit_format.lower():

            num_of_tus = ''
            greater_than = False
            less_than = False

            for char in number_of_tus:
                if char == '<':
                    less_than = True
                elif char == '>':
                    greater_than = True
                elif char.isdigit():
                    num_of_tus += char

            seq = []
            for i in range(int(num_of_tus)):
                if less_than:
                    seq.append("zero-or-one(promoter_any then cds)")
                elif greater_than:
                    seq.append("one-or-one(promoter_any then cds)")
                else:
                    seq.append("promoter_any then cds")

            seq = ' then '.join(seq)
            seq = f"({seq})"
            goldbar['N'].append(seq)


        else:

            num_of_tus = ''
            greater_than = False
            less_than = False

            for char in number_of_tus:
                if char == '<':
                    less_than = True
                elif char == '>':
                    greater_than = True
                elif char.isdigit():
                    num_of_tus += char

            seq = []
            for i in range(int(num_of_tus)):
                if less_than:
                    seq.append("zero-or-one(promoter_any then zero-or-one(promoter_any) then one-or-more(rbs then cds) then terminator_any)")
                elif greater_than:
                    seq.append("one-or-more(promoter_any then zero-or-one(promoter_any) then one-or-more(rbs then cds) then terminator_any)")
                else:
                    seq.append("promoter_any then zero-or-one(promoter_any) then one-or-more(rbs then cds) then terminator_any")

            seq = ' then '.join(seq)
            seq = f"({seq})"
            goldbar['N'].append(seq)

    if L:
        # Add to goldbar
        goldbar['L'].append("(zero-or-more(any_except_terminator_leaky) then zero-or-one(terminator_any))")

        # add any_except_terminator_notleaky to categories
        categories["any_except_terminator_leaky"] = {"promoter": parts["promoter"],
                                                        "ribosomeBindingSite": parts["ribosomeBindingSite"],
                                                        "cds": parts["cds"],
                                                        "terminator": parts["terminator_notleaky"]}

    if P:
        # Add to goldbar
        goldbar['P'].append("(zero-or-more(promoter_any then zero-or-one(promoter_notroadblocking) then one-or-more(rbs then cds) then terminator_any))")

        # Add promoter_notroadblocking to categories
        categories["promoter_notroadblocking"] = {"promoter": parts["promoter_notroadblocking"]}

    if O:
        # Extract orthogonal dictionary
        orthogonal_library_dict = orthogonal_library.to_dict('list')


        for part in list(orthogonal_library_dict.keys()):

            # Add part to categories
            categories = add_part(part, parts, categories)

            # Add any_except_part to categories
            categories = any_except_part(part, parts, categories, isConcrete)

            for part1 in orthogonal_library[part].dropna().to_list():
                # Add to goldbar
                goldbar['O'].append(f"(zero-or-more(any_except_{part}and{part1}) then zero-or-one(({part} then zero-or-more(any_except_{part1})) or ({part1} then zero-or-more(any_except_{part}))))")

                # Add part1 to categories
                categories = add_part(part1, parts, categories)
                
                # Add any_except_part1 to categories
                categories = any_except_part(part1, parts, categories, isConcrete)

                # Add any_except_partandpart1 to categories
                categories = any_except_part1andpart2(part, part1, parts, categories, isConcrete)

    if I:
        # Extract pji dictionary
        pji_library_dict = pji_library.to_dict('list')


        for part in list(pji_library_dict.keys()):

            # Multiple Parts in header?
            if '+' in part:
                part = part.split('+')
                name1 = "and".join(part)
                categories = any_except_multiple_parts(part, parts, categories, isConcrete)
                other_parts = pji_library["+".join(part)].dropna().to_list()
            else:
                name1 = part
                categories = any_except_part(part, parts, categories, isConcrete)
                other_parts = pji_library[part].dropna().to_list()
                part = [part]

            # Add part to categories
            categories = add_part(part, parts, categories)

            # Add Other parts to categories
            if len(other_parts) == 1:
                categories = add_part(other_parts[0], parts, categories)
                name2 = other_parts[0]
            else:
                categories = add_part(other_parts, parts, categories)
                name2 = "and".join(other_parts)

            # add any_except_partsandother_parts to categoies
            name3 = "any_except_" + name1 + "and" + name2
            categories = any_except_multiple_parts(part+other_parts, parts, categories, isConcrete)

            # Add to goldbar
            goldbar['I'].append(f"(zero-or-more(zero-or-more(any_except_{name1} or (one-or-more({name1}) then {name3})) then zero-or-more({name1})))")
                
    if R:
        for part in parts['do not repeat']:
            # Add to goldbar
            goldbar['R'].append(f"(zero-or-more(any_except_{part}) then zero-or-one({part} then zero-or-more(any_except_{part})))")

            # Add part to categories
            categories = add_part(part, parts, categories)

            # Add any_except_part to categories
            categories = any_except_part(part, parts, categories, isConcrete)
            
    if T:
        # Extract pji dictionary
        together_library_dict = together_library.to_dict('list')

        for part in list(together_library_dict.keys()):

            # Multiple Parts in header?
            if '+' in part:
                part = part.split('+')
                name1 = "and".join(part)
                other_parts = together_library["+".join(part)].dropna().to_list()
            else:
                name1 = part
                other_parts = together_library[part].dropna().to_list()
                part = [part]

            # Add part to categories
            categories = add_part(part, parts, categories)

            # Add Other parts to categories
            if len(other_parts) == 1:
                categories = add_part(other_parts[0], parts, categories)
                name2 = other_parts[0]
            else:
                categories = add_part(other_parts, parts, categories)
                name2 = "and".join(other_parts)

            # add any_part_except_parts to categoies
            name3 = "any_except_" + name1 + "and" + name2
            categories = any_except_multiple_parts(part+other_parts, parts, categories, isConcrete)


            # Add to goldbar
            seq = f"((zero-or-more({any_part_format}) then {name1} then zero-or-more({any_part_format}) then {name2} then zero-or-more({any_part_format}))"         
            seq += f" or (zero-or-more({any_part_format}) then {name2} then zero-or-more({any_part_format}) then {name1} then zero-or-more({any_part_format}))"
            seq += f" or (zero-or-more({name3})))"
            goldbar['T'].append(seq)
   
    if M:
        for part in parts['must include']:
            # Add to goldbar
            goldbar['M'].append(f"(zero-or-more({any_part_format}) then {part} then zero-or-more({any_part_format}))")

            # Add part to categories
            categories = add_part(part, parts, categories)

    if E:
        # Get Part
        part = parts['end_part']

        # Add part to categories
        categories = add_part(part, parts, categories)

        # Add any_except_part
        any_except_part(part, parts, categories, isConcrete)

        # Add to goldbar
        goldbar['E'].append(f"(zero-or-more(any_except_{part}) then one-or-more({part}))")

    # Join goldbar together with "and0"
    if isConcrete:
        and_tolerance = 'and0'
    else:
        and_tolerance = 'and1'

    # Make List of GOLDBAR
    goldbar_list = []
    for key in goldbar:
        for g in goldbar[key]:
            goldbar_list.append(g)
    
    goldbar_sequence = f" {and_tolerance} ".join(goldbar_list)

    # print goldbar and categories
    print(f'\nGoldbar:\n{goldbar_sequence}')
    print(f'\nCategories:')
    categories_json = json.dumps(categories)
    print(categories_json)

    write_txt = input('\nOutput Goldbar and Categories as txt file? (y/n): ')

    if write_txt.lower() == 'y':
        filename = input('Enter name for file (without .txt): ') + '.txt'
        print(f'Filename: {filename}')
        f = open(filename, "w")

        f.write("Goldbar:\n")
        f.write(goldbar_sequence)
        f.write('\n\n')

        f.write("Goldbar (rule by rule):\n")
        for i, key in enumerate(goldbar):
            if len(goldbar[key]) != 0:
                f.write(f"{key}: \n")
            for j, g in enumerate(goldbar[key]):
                f.write(g)
                if (i != len(goldbar)-1) and (j != len(goldbar[key]) -1):
                    f.write(f" {and_tolerance}\n")
                else:
                    f.write('\n')

        f.write("\n\n")
        categories_split = categories_json.split('},')
        f.write("Constellation Categories:\n")
        for i in range(len(categories_split)):
            f.write(categories_split[i])
            if i != len(categories_split) - 1:
                f.write('},\n')
        f.close()
        print('\nFile Created!')
        print('\nProgram Finished!')
    else:
        print('\nProgram Finished!')


if __name__ == "__main__":
    FILEPATH = input('Enter path to Part Library: ')
    part_library = pd.read_csv(FILEPATH, sep=',')
    print(part_library)

    circuit_format = input('\nEnter the format to use for genetic circuits: '
                        '\n"cello" - cello format (promoters, CDS, and terminators)'
                        '\n"PC" - promoters and CDS/Cassettes'
                        '\nNothing for default format (promoters, RBS, CDS, and terminators)'
                        '\nEnter: ')

    any_part_format = input('\n Enter whether you want parts to be "Concrete" or "Abstract":'
                            '\n"Concrete" - Parts with specific component IDs'
                            '\n"Abstract" - Parts without specific component IDs'
                            '\nDefault is Abstract'
                            '\nEnter: ')
    
    allIndividual = input('\nEnter if you want all individual parts in categories: '
                          '\n"n" - No, do not put all individual parts'
                          '\n"y" - Yes, put all individal parts'
                          '\nDefault is Yes'
                          '\nEnter: ')

    number_of_tus = input('\nEnter number of transcriptional units for each sequence: ')

    principles = input('\nEnter principles to use separated by a space [P L O I R]: '
                    '\n(P) Promoter Road Blocking: Road-blocking promoters must be alone or infront of regular promoters '
                    '\n(L) Leaky Terminators: Leaky Terminators can only be at end of circuit '
                    '\n(O) Not Orthogonal: Certain parts can not be in the same circuit '
                    '\n(I) Part-junction Interference: Certain parts can not be next to each other '
                    '\n(R) Do not Reuse Certain Parts: Certain Parts will not be reused in the same sequence'
                    '\n(T) Certain Parts Must be in the same circuit (header + 1 part in column)'
                    '\n(M) Parts that have to be included'
                    '\nEnter: ')

    goldbar_generator(principles, number_of_tus, part_library, circuit_format, any_part_format, allIndividual)
