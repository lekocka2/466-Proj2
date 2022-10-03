#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 15:44:32 2021

@author: logankocka
"""

#get the arguments from the command line call, save in list to reference later
import sys
import re
import matplotlib.pyplot as plt
# import numpy as np
    
# ------------ MODULES -------------------
def readFASTA(fileName, seqChoice):
    f = open(fileName, 'r')
    print(type(f))
    # allLines = f.readlines()
    # name is just accession number, description is everything else
    
    sequence = ""
    seq_name_list = []
    sequence_list =[]
    description_list = []
    
    if "fasta" in fileName:
        for line in f:
            if line[0] == ">":
                
                seq_name = line.split(" ",1)[0]
                seq_name = seq_name.strip('>')
                seq_name_list.append(seq_name)
                description = line.split(" ",1)[1]
                description = description.replace('\n','')
                description_list.append(description)
                if sequence:
                    sequence = sequence.replace("\n","")
                    sequence_list.append(sequence)
                    sequence = ''
            else:
                sequence += line
        
        sequence_list.append(sequence.replace("\n",""))
        sequence = sequence_list[seqChoice-1]
        seq_name = seq_name_list[seqChoice-1]
        description = description_list[seqChoice-1]
        
    # if it is a fastq file run this 
    if "fastq" in fileName:
        allLines = open(fileName, 'r').readlines()
        for i in range(0,len(allLines),4):
            seq = allLines[i:i+4]
            l = seq[:2]
            seq_name = l[0].split(" ")[0]
            description = l[0].split(" ")[1]
            sequence = l[1].replace("\n","")
   
    
    return seq_name,description,sequence,sequence


def printlnFASTA(seqName, seqDescription, sequence, n):
    # print the name and description
    print(">" + seqName + " " + seqDescription)
    # FASTA_LINE_NUM = n
    # for i in range(0,len(sequence),FASTA_LINE_NUM):
    #     print(sequence[i:i+FASTA_LINE_NUM])
    for i in range(0, len(sequence),n):
        print(sequence[i:i+n])
    return


def printWithRuler(seqName, sequence, seqDescription, spacerChoice, n, indexes):
    if indexes == "":
        NON_FASTA_LINE_NUM = n
        if spacerChoice == "Y":
            spacer = " "
        if spacerChoice == "N":
            spacer = ""
        
        # print first line and numbers with spaces
        print('>' + seqName + ' ' + seqDescription + '\n')
        print(15*' ', end="")
        nums = list(range(1,11))
        for i in nums:
            print(nums[i-1], end="")
            if i < 9:
                print(9*' ' + spacer, end="")
            if i == 9:
                print(8*' ' + spacer, end="")
        print('\n' + " Line", end=" ")
        for i in nums:
            print("1234567890" + spacer, end="")
            if i == 10: print()
        
        # print each line with spacer every ten characters
        if spacerChoice == "N":
            # format the list
            NON_FASTA_LINE_NUM = 100 # split every 100 characters
            l_N = [sequence[i:i+NON_FASTA_LINE_NUM] for i in range(1, len(sequence), NON_FASTA_LINE_NUM)]
            #print them
            for count,val in enumerate(l_N):
                print('{:>5} {:<100}'.format(count+1,val))
        
        if spacerChoice == "Y":
            # format the list
            l_Y = [sequence[i:i+10] for i in range(1, len(sequence), 10)]
            l_Y2 = []; start = 0
            for i in l_Y:
                str = " ".join(l_Y[start:start+10])
                l_Y2.append(str)
                start += 10
            l_Y3 = list(filter(None, l_Y2))
            # print them
            for count,val in enumerate(l_Y3):
                print('{:>5} {:<110}'.format(count+1,val))

    else:
        NON_FASTA_LINE_NUM = n
        if spacerChoice == "Y":
            spacer = " "
        if spacerChoice == "N":
            spacer = ""
        
        # print first line and numbers with spaces
        print('>' + seqName + ' ' + seqDescription + '\n')
        print(15*' ', end="")
        nums = list(range(1,11))
        for i in nums:
            print(nums[i-1], end="")
            if i < 9:
                print(9*' ' + spacer, end="")
            if i == 9:
                print(8*' ' + spacer, end="")
        print('\n' + " Line", end=" ")
        for i in nums:
            print("1234567890" + spacer, end="")
            if i == 10: print()
        
        # find and replace the motifs with lower case letters
        seq_list = list(sequence)
        for i in indexes:
            # for j in range(indexes[i][0],indexes[i][-1]):
            
            seq_list[i[0]:i[1]] = [j.lower() for j in seq_list[i[0]:i[1]]]
            
        print("TEST", ''.join(seq_list))
        sequence = ''.join(seq_list)
        # print each line with spacer every ten characters
        if spacerChoice == "N":
            # format the list
            NON_FASTA_LINE_NUM = 100 # split every 100 characters
            l_N = [sequence[i:i+NON_FASTA_LINE_NUM] for i in range(1, len(sequence), NON_FASTA_LINE_NUM)]
            #print them
            for count,val in enumerate(l_N):
                print('{:>5} {:<100}'.format(count+1,val))
        
        if spacerChoice == "Y":
            # format the list
            l_Y = [sequence[i:i+10] for i in range(1, len(sequence), 10)]
            l_Y2 = []; start = 0
            for i in l_Y:
                str = " ".join(l_Y[start:start+10])
                l_Y2.append(str)
                start += 10
            l_Y3 = list(filter(None, l_Y2))
            # print them
            for count,val in enumerate(l_Y3):
                print('{:>5} {:<110}'.format(count+1,val))


    return


def nucleotideCounter(sequence):
    # initialize
    counts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
  
    for i in sequence: # count values
        if i in counts:
            counts[i] += 1
        else:
            counts[i] = 1
    
    print("(2.1) Nucleotide Counts: ") # print the key value pairs
    if 'N' in counts: print("A=[{A}] T=[{T}] G=[{G}] C=[{C}] N=[{N}]".format(**counts))
    else: print("A=[{A}] T=[{T}] G=[{G}] C=[{C}]".format(**counts))
    # some don't have N, some do...
    return counts


def gcContent(sequence):
    # find percentage of all nucleotides are G or C
    counts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
    for i in sequence: # count values
        if i in counts:
            counts[i] += 1
        else:
            counts[i] = 1
    total = sum(counts.values())
    
    GC_Counts = int(counts.get('G')) + int(counts.get('C'))
    
    return print("(2.2) GC content: ", round(GC_Counts/total*100,2), "%")


def diCounter(sequence):
    # declare an empty dictionary
    di_dict = {}
    # initiate with 16 dinucleotides as the keys with values of 0
    ### non hard coded, use loops to initiate the slices and make counts zero
    
    # for i in range(len(sequence)):
    #     # if i < len(sequence):
    #     slice = sequence[i:i+2:1]
    #     if slice not in di_dict:
    #         if len(slice) == 2:
    #             di_dict[slice] = 0
                
    let = ['T', 'C', 'A', 'G']
    newList = []
    for l in let:
        for l2 in let:
            newList += [l + l2]
    for i in newList:
        di_dict[i] = 0
    
    for i in range(len(sequence)):
        if i != (len(sequence)):
            slice = sequence[i:i+2:1]
            if len(slice) == 2:
                if not 'N' in slice:
                    di_dict[slice] += 1
                    
    # print dynamically
    print("(2.3) Di-nucleotide Counts: ")
    print("TT={:<3} TC={:<3} TA={:<3} TG={:<3}".format(di_dict['TT'],di_dict['TC'],di_dict['TA'],di_dict['TG']))
    print("CT={:<3} CC={:<3} CA={:<3} CG={:<3}".format(di_dict['CT'],di_dict['CC'],di_dict['CA'],di_dict['CG']))
    print("AT={:<3} AC={:<3} AA={:<3} AG={:<3}".format(di_dict['AT'],di_dict['AC'],di_dict['AA'],di_dict['AG']))
    print("GT={:<3} GC={:<3} GA={:<3} GG={:<3}".format(di_dict['GT'],di_dict['GC'],di_dict['GA'],di_dict['GG']))
    
    return di_dict


def codonProfile(sequence):
    # takes in DNA sequence, returns dictionary of 64 keys and zero values
    codon_dict = {} #declare empty
    
    let = ['T', 'C', 'A', 'G']
    newList = []
    for l in let:
        for l2 in let:
            for l3 in let:
                newList += [l + l2 + l3] # get combos of all letters
    #add these to dictionary
    for i in newList:
        codon_dict[i] = 0 # initialize all with zero as value

    return codon_dict


def printCodonProfile(sequence):
    dict = codonProfile(sequence)
    # lets = ['T', 'C', 'A', 'G']
    codon_order = ['TTT', 'TCT', 'TAT', 'TGT', 'TTC', 'TCC', 'TAC', 'TGC',
                   'TTA', 'TCA', 'TAA', 'TGA', 'TTG', 'TCG', 'TAG', 'TGG',
                   'CTT', 'CCT', 'CAT', 'CGT', 'CTC', 'CCC', 'CAC', 'CGC',
                   'CTA', 'CCA', 'CAA', 'CGA', 'CTG', 'CCG', 'CAG', 'CGG',
                   'ATT', 'ACT', 'AAT', 'AGT', 'ATC', 'ACC', 'AAC', 'AGC',
                   'ATA', 'ACA', 'AAA', 'AGA', 'ATG', 'ACG', 'AAG', 'AGG',
                   'GTT', 'GCT', 'GAT', 'GGT', 'GTC', 'GCC', 'GAC', 'GGC',
                   'GTA', 'GCA', 'GAA', 'GGA', 'GTG', 'GCG', 'GAG', 'GGG']
    
    for i in range(len(sequence)):
        if i != (len(sequence)):
            slice = sequence[i:i+3:1]
            if len(slice) == 3:
                if not 'N' in slice:
                    dict[slice] += 1
    
    print('\n\n(2.4) Codon Profile')
    
    print("\t\t      2nd")
    print("-"*45)
    print("{0}    {1}       {2}       {3}       {4}          {5}".format('1st', 'T','C','A','G', '3rd'))
    for i in range(0,64,4):
        if i==0:
            print('T     ', end="")
        elif i==16:
            print('C     ', end="")
        elif i==32:
            print('A     ', end="")
        elif i==48:
            print('G     ', end="")
        else: print('      ', end="")
        
        for j in [0,1,2,3]:
            my_codon = codon_order[i+j]
            print("{0}={1}".format(my_codon, repr(dict[my_codon]).rjust(3)), end=' ')
            # if codon_order.index(j)+1%8 == 0:
            if j == 3:
                print('   ', codon_order[i+j][-1], end="\n")
    print('\n')
    return


def codonProfileCompare(seq1, seq2):
    dict1 = codonProfile(seq1)
    dict2 = codonProfile(seq2)
    
    codon_order = ['TTT', 'TCT', 'TAT', 'TGT', 'TTC', 'TCC', 'TAC', 'TGC',
                   'TTA', 'TCA', 'TAA', 'TGA', 'TTG', 'TCG', 'TAG', 'TGG',
                   'CTT', 'CCT', 'CAT', 'CGT', 'CTC', 'CCC', 'CAC', 'CGC',
                   'CTA', 'CCA', 'CAA', 'CGA', 'CTG', 'CCG', 'CAG', 'CGG',
                   'ATT', 'ACT', 'AAT', 'AGT', 'ATC', 'ACC', 'AAC', 'AGC',
                   'ATA', 'ACA', 'AAA', 'AGA', 'ATG', 'ACG', 'AAG', 'AGG',
                   'GTT', 'GCT', 'GAT', 'GGT', 'GTC', 'GCC', 'GAC', 'GGC',
                   'GTA', 'GCA', 'GAA', 'GGA', 'GTG', 'GCG', 'GAG', 'GGG']
    
    for i in range(len(seq1)):
        if i != (len(seq1)):
            slice1 = seq1[i:i+3:1]
            if len(slice1) == 3:
                if not 'N' in slice1:
                    dict1[slice1] += 1
                    
    for i in range(len(seq2)):
        if i != (len(seq2)):
            slice2 = seq2[i:i+3:1]
            if len(slice2) == 3:
                if not 'N' in slice2:
                    dict2[slice2] += 1
    
    #create new dictionary comparing the two for printing
    comp_dict = {}
    for i in codon_order:
        if dict1[i] > dict2[i]:
            if dict1[i] - dict2[i] < 10:
                comp_dict[i] = "   "
            elif dict1[i] - dict2[i] < 20:
                comp_dict[i] = " > "
            elif dict1[i] - dict2[i] < 30:
                comp_dict[i] = " >>"
            else:
                comp_dict[i] = ">>>"
        else:
            if dict2[i] - dict1[i] < 10:
                comp_dict[i] = "   "
            elif dict2[i] - dict1[i] < 20:
                comp_dict[i] = " > "
            elif dict2[i] - dict1[i] < 30:
                comp_dict[i] = " >>"
            else:
                comp_dict[i] = ">>>"
    
    print("\t\t      2nd")
    print("-"*45)
    print("{0}    {1}       {2}       {3}       {4}          {5}".format('1st', 'T','C','A','G', '3rd'))
    for i in range(0,64,4):
        if i==0:
            print('T     ', end="")
        elif i==16:
            print('C     ', end="")
        elif i==32:
            print('A     ', end="")
        elif i==48:
            print('G     ', end="")
        else: print('      ', end="")
        
        for j in [0,1,2,3]:
            my_codon = codon_order[i+j]
            print("{0}={1}".format(my_codon, comp_dict[my_codon].rjust(3)), end=' ')
            if j == 3:
                print('   ', codon_order[i+j][-1], end="\n")
    
    # print the plots here
    # use matplotlib  
    plt.subplot(2,1,1)
    plt.bar(dict1.keys(),dict1.values())
    plt.ylabel('Frequency')
    plt.xticks([],[])
    plt.title('Seq 1')
    plt.subplot(2,1,2)
    plt.bar(dict2.keys(),dict2.values())
    plt.ylabel('Frequency')
    plt.xticks(rotation='vertical')
    plt.title('Seq 2')
    plt.show()
    
    return


def inquiry(sequence, frag):
    start = int(frag.partition('::')[0])
    end = int(frag.partition('::')[2])
    
    #print fragment length
    print("The fragment you selected has a length of ", end-start, " nucleotides:")
    numDashes = (end-start) - 1 - (int(len(str(start)))) - (int(len(str(end))))
    print("<",start,"-"*numDashes,end,">", sep="")
    seqFrag = sequence[start:end+1]
    print("|"*len(seqFrag))
    print(seqFrag) # keep this
    nucleotideCounter(seqFrag)
    if 'G' in seqFrag and 'C' in seqFrag:
        gcContent(seqFrag)
    elif 'G' in seqFrag or 'C' in seqFrag:
        print("GC content: ", round((1/len(seqFrag)),2), "%")
    else:
        print("GC content: 0%")
    
    diCounter(seqFrag)
    
    return 


def translation(sequence):
    sequence = sequence.replace('N', "") # remove Ns
    rnaseq = sequence.replace('T', 'U')
    prot_pairs = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 
                  'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
                  'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', 
                  'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 
                  'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 
                  'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 
                  'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 
                  'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 
                  'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*', 
                  'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 
                  'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 
                  'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 
                  'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W', 
                  'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 
                  'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 
                  'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
                  }
    prot_seq = ""
    for i in range(0,len(rnaseq),3):
        slice = rnaseq[i:i+3]
        if len(slice) == 3:
            prot_seq += prot_pairs[slice]
    return prot_seq


def ORF(sequence):
    sequence = sequence.replace('N', "") # remove Ns
    
    base_pairs = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    rev_seq_ = ""
    for i in sequence:
        rev_seq_ += base_pairs[i]
    
    # shift and translate
    ps1 = translation(sequence)
    ps2 = translation(sequence[1:])
    ps3 = translation(sequence[2:])
    ps_1 = translation(rev_seq_)
    ps_2 = translation(rev_seq_[:len(rev_seq_)-1])
    ps_3 = translation(rev_seq_[:len(rev_seq_)-2])
    # find motifs using regex
    frame1 = re.findall('M\w{10,}\*', ps1)
    frame2 = re.findall('M\w{10,}\*', ps2)
    frame3 = re.findall('M\w{10,}\*', ps3)
    frame_1 = re.findall('M\w{10,}\*', ps_1)
    frame_2 = re.findall('M\w{10,}\*', ps_2)
    frame_3 = re.findall('M\w{10,}\*', ps_3)
    
    empty = "##### NO VALID ORF #####"
    print("(2.5) 6-frame translations: \n") 
    # frame +1
    if not frame1: 
        print("ORF (+1): ", empty)
    elif len(frame1) == 1: 
        print("ORF (+1): ", "".join(frame1))
    else: 
        for i in sorted(frame1,key=len):
            print("ORF (+1): ", i)
    # frame +2
    if not frame2: 
        print("ORF (+2): ", empty)
    elif len(frame2) == 1: 
        print("ORF (+2): ", "".join(frame2))
    else: 
        for i in sorted(frame2,key=len):
            print("ORF (+2): ", i)
    # frame +3 
    if not frame3: 
        print("ORF (+3): ", empty)
    elif len(frame3) == 1: 
        print("ORF (+3): ", "".join(frame3))
    else: 
        for i in sorted(frame3,key=len):
            print("ORF (+3): ", i)
    # frame -1
    if not frame_1: 
        print("ORF (-1): ", empty)
    elif len(frame_1) == 1: 
        print("ORF (-1): ", "".join(frame_1))
    else: 
        for i in sorted(frame_1,key=len):
            print("ORF (-1): ", i)
    # frame -2
    if not frame_2: 
        print("ORF (-2): ", empty)
    elif len(frame_2) == 1: 
        print("ORF (-2): ", "".join(frame_2))
    else: 
        for i in sorted(frame_2,key=len):
            print("ORF (-2): ", i)
    # frame -3
    if not frame_3: 
        print("ORF (-3): ", empty)
    elif len(frame_3) == 1: 
        print("ORF (-3): ", "".join(frame_3))
    else: 
        for i in sorted(frame_3,key=len):
            print("ORF (-3): ", i)
    print()
    return 


def motif_finder(motif_list, sequence):     
    dict = {} # initalize empty dict
    for i in motif_list: # for i in the input list
        # print(i)
        result = re.findall(i,sequence)
        
        dict[i] = [len(result)] # dict with key motif = length of the result list/num matches
        p = re.compile(i) 
        iter = p.finditer(sequence) 
        for match in iter:
            dict[i].append(match.span())
    # print(dict)
    
    return dict

# ------- MAIN --------------------------------
def main():
    args = sys.argv
    fileName = args[1]
    user = args[0].split("_")[0]
    
    #print the welcome message
    print("\nWelcome Sequence Viewer! Programmer: ", user)
    #print sequences detected
    num = len([1 for line in open(fileName) if line.startswith(">")])
    if num == 1:
        print("There is 1 sequence detected in the file: ", fileName)
    else: 
        print("There are ", num, " sequences detected in the file: ", fileName)

    print("Which sequence do you want to examine ", [*range(1, num+1)], "?")
    seqChoice = int(input())

    print("\n[Part I]: Display Mode \n")
    
    # get viewing format choice
    while True:
        try:
            formatChoice = input("Do you want to view the sequence in FASTA format? (Y/N) ")
            if formatChoice == "Y" or formatChoice == "N":
                break;
            else: 
                print("Invalid. Enter Y or N.")
        except:
            continue

    # read the file
    new_lst = readFASTA(fileName, seqChoice)

    # if yes, print Screen 2
    if formatChoice == "Y":
        printlnFASTA(new_lst[0], new_lst[1], new_lst[2], 60)

    #if N, ask about spacer option
    if formatChoice == "N":
        while True:
            try:
                spacerChoice = input("Do you need a spacer for viewing nucleotide positions? (Y/N) ")
                if spacerChoice == "Y" or spacerChoice == "N":
                    break;
                else: 
                    print("Invalid. Enter Y or N.")
            except:
                continue
    
    if formatChoice == "N":
        # print with or without spacer
        printWithRuler(new_lst[0], new_lst[2], new_lst[1], spacerChoice, 100, 'N')
        
    print("\n[Part II]: Analysis Mode \n")
    nucleotideCounter(new_lst[2])
    gcContent(new_lst[2])
    diCounter(new_lst[2])
    printCodonProfile(new_lst[2])
    ORF(new_lst[2]) 
    
    print("[Part III]: Inquiry Mode \n")
    print("(3.1) Extract a DNA fragment \n(3.2) Compare two codon profiles \n(3.3) Find a motif\n")
    option = input("Your choice is: ")
    
    #loop for getting the right kind of input for viewing a fragment
    while True:
        if option == "3.1":
            try:
                frag = input("Please enter the start and end positions (e.g., 19::48): ")
                # if frag == "N":
                #     sys.exit()
                if "::" not in frag:
                    print("Separate start and end value with ""::""")
                    continue
                start = int(frag.partition('::')[0])
                end = int(frag.partition('::')[2])
                # if a fragment is entered, call the function 
            except ValueError:
                print("Start and end values must be integers.")
                continue
            if start < 0:
                print("Start value must be nonnegative.")
                continue
            elif end > len(new_lst[2]):
                print("End value can't be larger than sequence length.")
                continue
            elif start-end > 0:
                print("End value must be higher than start value.")
                continue
            else:
                inquiry(new_lst[2], frag)
                another = input("Do you want to extract another fragment? ")
                if another == "N":
                    sys.exit()
                continue
        
        elif option == "3.2":
            # file must be in same directory as this script to work
            try: 
                secFile = input("Please open another sequence file that contains only 1 sequence.\nEnter the file name here: ")
                break 
                # if not os.path.exists(secFile):
                #     print("That file does not exist in the current directory.")
            except Exception as e:
                print(e)
            finally:
                # print the two file names
                print("File 1: ", fileName)
                print("File 2: ", secFile)
                # call the function
                secFileSeq = readFASTA(secFile, 1)[2]
                codonProfileCompare(new_lst[2],secFileSeq)
                # when all done, exit
                sys.exit() 
                    
                
        elif option == "3.3":
            try:
                #get the motif seq
                print("Please enter the motif that you want to find (e.g., AATTC,GCCTTA)")
                motif_list = [i.strip() for i in input("Enter the motif sequence here: ").split(',')]
                break
            except Exception as e:
                print(e)
            finally:
                # implement motif_finder function
                dict = motif_finder(motif_list, new_lst[2])
                # this returns a dictionary with all info needed
                # num = len(max(motif_list))
                print("{:<15}  {:<8}  {:<19}".format("Motif", "Frequency", "Position(Start-End)"))
                for i,k in dict.items():
                    # print the motif, freq, and indexes in each row
                    print("{:<15}  {:<11}".format(i, k[0]), end="")
                    for i in range(1,k[0]+1):
                        print("{:<11}".format(repr(k[i])),end="")
                    print()
                print()
                
                printWithRuler(new_lst[0], new_lst[2], new_lst[1], 'Y', 100, k[1:])
                
                # when all done, exit
                sys.exit()
    
    return


if __name__=='__main__':
    main()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    