# David Bogoslavsky
# Date: 9/1/22
import csv
import os
import sys
import random
import re
import json


class Polymerase:
    """
    Polymerase is responsible for transcribe a dna_seq into the complementary rna/dna seq.
    """
    def __init__(self, type, error_rate=0.0):
        """
        constructor of Polymerase
        :param type: string for the type of pol
        :param error_rate:  float of the number of errors.
        """
        self.type = type
        self.error_rate = error_rate
        self.num_of_errors = 0

    def generating_errors(self, type_seq, num_of_errors, genetic_seq_string):
        """
        This method responsible for generating errors of polymerase
        :param type_seq: if DNA of RNA
        :param num_of_errors: how
        :param genetic_seq_string: the genetic code for polymerase to work on.
        :return: genetic_seq_string_returned - The seq out from the polymerase that can be RNA, or DNA.
        """
        if type_seq == 'RNA':
            list_of_nucleotides = ['A', 'U', 'G', 'C']
        elif type_seq == 'DNA':
            list_of_nucleotides = ['A', 'T', 'G', 'C']
        else:
            raise AssertionError(type_seq)
        list_of_already_replaced_nucs = []
        genetic_seq = []
        index_of_error_nucleotide = True
        for letter in genetic_seq_string:
            genetic_seq.append(letter)  # converting string to list for easier work
        while num_of_errors != 0:  # we generate the number of errors as calculated via error_rate before
            index_was_already_replaced = True
            while index_was_already_replaced:  # we check whether we already generated an error in that nuc
                index_of_error_nucleotide = random.randint(0, len(genetic_seq) - 1)
                if list_of_already_replaced_nucs:
                    for index in list_of_already_replaced_nucs:
                        if index == index_of_error_nucleotide:
                            index_was_already_replaced = True  # if we already generated an error, we will call diff nuc
                            break
                        else:
                            index_was_already_replaced = False
                else:
                    index_was_already_replaced = False
                if index_was_already_replaced:
                    continue
                list_of_already_replaced_nucs.append(index_of_error_nucleotide)  # we mark for next this
                # nuc as errored
            errored_nuc = genetic_seq[index_of_error_nucleotide]
            list_of_nucleotides.remove(errored_nuc)  # we need to replace the nuc with other nuc for it to be errored
            replacing_errored_nuc = random.randint(0, 2)
            errored_nuc_replaced = list_of_nucleotides[replacing_errored_nuc]  # picking replacing nuc
            list_of_nucleotides.append(errored_nuc)
            genetic_seq[index_of_error_nucleotide] = errored_nuc_replaced  # creating the error in our seq!
            num_of_errors -= 1  # marking that we've created one error
        genetic_seq_string_returned = ''.join(genetic_seq)  # converting back to string
        return genetic_seq_string_returned

    def transcribe(self, dna_seq):
        """
          DNA/RNA polymerase function.
          :param dna_seq: the input of seq
          :return: new seq.
          """
        if dna_seq:
            self.num_of_errors = len(dna_seq) * self.error_rate  # every time the number is restarted when the
            # function is called
            is_the_number_not_round = isinstance(self.num_of_errors, float)
            if is_the_number_not_round and self.num_of_errors != int(self.num_of_errors):  # we check if it's not
                # round, and if its not exactly 1.0
                self.num_of_errors += 1  # if the number is not round, we add 1,
                # so that the conversion to int will result in
                # rounding up!
                self.num_of_errors = int(self.num_of_errors)
            self.num_of_errors = int(self.num_of_errors)
            if self.type == 'RNA':
                dna = dna_seq.upper()
                rna = ''  # initialize RNA to empty string
                dna_to_rna_dict = {"A": "U", "T": "A", "C": "G", "G": "C"}
                for letter in dna:  # loop over DNA sequence characters to make them complementary nucleotide
                    if letter in dna_to_rna_dict:
                        rna += dna_to_rna_dict[letter]
                    else:
                        raise AssertionError(dna_seq)
                        # rna = None
                        # break
                if rna:
                    rna = rna[::-1]  # reverse the string
                    if self.num_of_errors > 0:
                        rna = self.generating_errors('RNA', self.num_of_errors, rna)
                return rna  # return output RNA sequence
            elif self.type == 'DNA':
                dna = dna_seq.upper()
                rep_dna = ''  # initialize dna to empty string
                dna_to_dna_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
                for letter in dna:  # loop over DNA sequence characters to make them complementary nucleotide
                    if letter in dna_to_dna_dict:
                        rep_dna += dna_to_dna_dict[letter]
                    else:
                        raise AssertionError(dna_seq)
                        # rep_dna = None
                        # break
                if rep_dna:
                    rep_dna = rep_dna[::-1]  # reverse the string
                    if self.num_of_errors > 0:
                        rep_dna = self.generating_errors('DNA', self.num_of_errors, rep_dna)
                return rep_dna  # return output DNA sequence
            else:
                ret_string = ''
                return ret_string
        else:
            raise AssertionError(dna_seq)
            # return None


class Ribosome:
    """
    Ribosome Class: the ribosome responsible on synthesizing rna seq into a protein.
    """
    def __init__(self, genetic_code, start_codons):
        """
        Initialize Ribosome - constructor
        :param genetic_code: a map of codons and the letter of protein
        :param start_codons: a set of codons that means start translation for ribosome.
        """
        self.genetic_code = genetic_code
        self.start_codons = start_codons
        # edit afterwards
        self.stop_codons = []
        for codon in self.genetic_code:
            if self.genetic_code[codon]:
                pass
            else:
                self.stop_codons.append(codon)
        # end of adding

    def translate(self, rna_seq):
        """
        A sub unit of synthesize function:
        able to translate the rna seq to a list of codons for the use of ribosome
        :return: list type - the longest reading frame
        """
        if rna_seq:
            rna = rna_seq.upper()
            ret = None
            rna_longest_seq = ''
            for start_codon in self.start_codons:
                for cur_shift in range(0, 3):  # we check all 3 possible frame shifts
                    remaining_sequence = rna
                    while len(remaining_sequence) > len(
                            rna_longest_seq):  # we do this loop in case we have an end codon in frame
                        # we continue checking for longer proteins after the current long seq
                        counter = 3  # we take the AUG in consideration
                        while True:
                            idx = remaining_sequence.find(start_codon)  # find start index of Methionine codon in rna
                            if idx == -1:  # if not found
                                break
                            if idx % 3 == cur_shift:  # it means that our found AUG is in our frame shift
                                break
                            remaining_sequence = remaining_sequence[idx + 3:]
                            # if the start codon we found is not in our frameshift we check the next remaining codons
                        if idx == -1:
                            # if we haven't found AUG, we break the loop, because what is next is not relevant
                            break
                        end_codon_found = False
                        if len(rna) == 3:  # a case for a single protein which is AUG only
                            rna_longest_seq = start_codon
                        for codon_start_index in range(idx, len(remaining_sequence) - 2, 3):
                            # we check and create a sequence of nucleotides
                            cur_codon = remaining_sequence[codon_start_index: codon_start_index + 3]
                            counter += 3
                            if counter > len(rna_longest_seq):
                                rna_longest_seq = remaining_sequence[idx: codon_start_index + 3]
                                # we change the biggest rna seq if we find that other frame shift has a longer one
                            # if cur_codon == 'UAG' or cur_codon == 'UAA' or cur_codon == 'UGA':  # end codons
                            if cur_codon in self.stop_codons:
                                if counter > len(rna_longest_seq):
                                    rna_longest_seq = remaining_sequence[idx: codon_start_index]
                                # we don't include the end codon
                                end_codon_found = True  # declare we found end codon in frame
                                # remaining_sequence = remaining_sequence[
                                #                     codon_start_index + 3 - cur_shift:]
                                remaining_sequence = remaining_sequence[
                                                     idx + 3:]  # check seq after end
                                break
                        if end_codon_found:
                            # if we found end codon, we continue looping in the same frame, after the end codon
                            pass
                        else:  # if else, we move to the next frame
                            break
            if len(rna_longest_seq) >= 3:
                sent_rna = []
                for i in range(0, len(rna_longest_seq), 3):
                    next_codon = rna_longest_seq[i:i + 3]
                    if len(next_codon) != 3:  # we cut remaining nucleotides if they are not 3
                        break
                    sent_rna.append(next_codon)
                ret = sent_rna
            return ret
        else:
            raise AssertionError(rna_seq)
            # return None

    def synthesize(self, rna_seq):
        """
        making list of amino-acids by genetic mRNA - rna_seq.
        :param rna_seq: the rna for translation and synthesize to a neckless of amino-acids
        :return: a string of protein from ribosome
        """
        ret_protein_seq = ''  # the seq we will returned - protein
        translated_seq = self.translate(rna_seq)
        protein_dict = self.genetic_code
        if translated_seq:  # if there is something to synthesise
            for codon in translated_seq:
                if codon in protein_dict:
                    # saving the proteins one after one in returned list
                    ret_protein_seq += protein_dict[codon]
                else:  # in case there is no more codons to make form them amino-acid:
                    break
        return ret_protein_seq


class Cell:
    """
    Cell class : general cell.
    """
    def __init__(self, name, genome, num_copies, genetic_code, start_codons, division_rate, error_rate=0.0):
        """
        Initialize new Cell
        :param name: string name (3 options)
        :param genome: list string of seqs.
        :param num_copies: natural number
        :param genetic_code: a one from two given maps
        :param start_codons: set of codons means start to ribosome
        :param division_rate: num between 0-1, default 0
        :param  error_rate: double, default 0
        """
        self.name = name  # string
        self.genome = genome  # list of dna seqs
        self.num_copies = num_copies  # an int number
        self.genetic_code = genetic_code  # a dict which keys are the codons and the values are the amino-acids
        self.start_codons = start_codons  # set
        self.division_rate = division_rate
        # cell's Ribosome:
        self.cell_ribosome = Ribosome(self.genetic_code, self.start_codons)
        self.error_rate = error_rate
        # cell's Polymerases:
        self.cell_DNA_polymerase = Polymerase('DNA', self.error_rate)
        self.cell_RNA_polymerase = Polymerase('RNA', self.error_rate)

    def __repr__(self):
        """
        for print cell in general
        :return: string value
        """
        return "<" + self.name + ", " + str(self.num_copies) + ", " + str(self.division_rate) + ">"

    def __mul__(self, num):
        """
        :param: num of cells to replicated
        :return: list of cells - by the other num as inserted
        """
        new_cell = Cell(self.name, self.genome, self.num_copies, self.genetic_code,
                        self.start_codons, self.division_rate)
        new_divided_cells = []
        for number_of_divided in range(num):  # num is the number of cells replicated
            new_divided_cells.append(new_cell)
        return new_divided_cells

    def mitosis(self):
        """
         :return: cell after mitosis.
         """
        mitosis_list = []
        for i in range(self.division_rate):
            cell_division = Cell(self.name, self.genome, self.num_copies, self.genetic_code, self.start_codons,
                                 self.division_rate)
            # add to new cells list:
            mitosis_list.append(cell_division)
        # in list - all new cells from mitosis.
        return mitosis_list

    def meiosis(self):
        """
        decrease division
        :return: two cells with half of the copies the older cell had
        """
        if self.num_copies % 2 == 0: #if even
            meiosis_list = [] #make meiosis list
            complementary_cell = [] #cell 
            for gene in self.genome: #for all genomes
                complementary_cell.append(self.cell_DNA_polymerase.transcribe(gene))
            meiosis_cell2 = Cell(self.name, complementary_cell, self.num_copies / 2,
                                 self.genetic_code, self.start_codons, self.division_rate)
            meiosis_cell1 = Cell(self.name, self.genome, self.num_copies / 2,
                                 self.genetic_code, self.start_codons, self.division_rate)
            meiosis_list.append(meiosis_cell1) #add cell1
            meiosis_list.append(meiosis_cell2) #add cell2
            return meiosis_list #return the list (2 cells - one complementy to other)
        else:
            return None

    def find_srr(self, dna_seq):
        """
            in 3 loops - go over all the possible options to reveal the srr
            :param dna_seq: seg of DNA
            :return: list of subs appears in their maximum appearance.
            """
        if dna_seq:
            dna_seq.upper()
            legal_seq = True
            ret_list = []
            seq_list = []
            counter_list = []
            for letter in dna_seq:  # loop over DNA sequence characters to check legality
                if letter == 'A':
                    continue
                elif letter == 'T':
                    continue
                elif letter == 'C':
                    continue
                elif letter == 'G':
                    continue
                else:
                    legal_seq = False
            if legal_seq:
                for short_seq_len in range(1, 7):  # we begin with sequences of 1 and end with 6
                    for start_index in range(0, len(dna_seq) - short_seq_len):
                        # we run through every sequence from the first letter to the last
                        short_seq = dna_seq[
                                    start_index:start_index + short_seq_len]  # the sequence we compare the next ones
                        if short_seq == dna_seq[start_index - short_seq_len:start_index]:
                            continue
                            # we checked above whether the cur seq equals to previous seq in distance of the same len
                        counter = 1
                        for i in range(start_index + short_seq_len, len(dna_seq) + 1 - short_seq_len, short_seq_len):
                            cur_seq = dna_seq[i:i + short_seq_len]  # the next sequences after short seq
                            if short_seq == cur_seq:  # we check whether the next sequence is equal to the previous
                                counter += 1  # we add 1 to the counter if so
                            if short_seq != cur_seq or i + short_seq_len * 2 > len(dna_seq):
                                # we've reached the final count. Either we had a mismatch or the string is over
                                if counter >= 3:
                                    # if our next sequence doesn't equal to before, we check 3 or more for printing.
                                    # the concept: we create three list with which we would work comfortably
                                    ret_list.append(short_seq + "," + str(counter))  # main list of srr
                                    seq_list.append(short_seq)  # list of the seqs without the "," and counter
                                    counter_list.append(counter)  # list of ordered counter numbers
                                    length = len(ret_list)
                                    for list_index in range(length - 1):
                                        if seq_list[length - 1] == seq_list[list_index]:
                                            # here we check for same seq,
                                            # and then check the higher counter and delete the redundant
                                            if counter_list[length - 1] > counter_list[list_index]:
                                                ret_list.remove(ret_list[list_index])
                                                seq_list.remove(seq_list[list_index])
                                                counter_list.remove(counter_list[list_index])
                                            else:
                                                ret_list.remove(ret_list[length - 1])
                                                seq_list.remove(seq_list[length - 1])
                                                counter_list.remove(counter_list[length - 1])
                                break
            ret_list.sort()
            ret_string = ''
            for srr in ret_list:
                ret_string += srr
                if srr != ret_list[len(ret_list) - 1]:
                    ret_string += ';'
            return ret_string
        else:
            raise AssertionError(dna_seq)
            # return None

    def repertoire(self):
        """
        for printing a repertoire for user - summery
        :return: tuple list
        """
        repertoire_list = []  # for saving the future tuple
        for gene in self.genome:  # on each genome
            srr_ret = self.find_srr(gene)  # find the srrs
            if srr_ret:  # if there is at list one srr in list
                pass
            else:
                srr_ret = 'No simple repeats in DNA sequence'
            rna_seq = Polymerase('RNA', 0).transcribe(gene)
            protein_seq = self.cell_ribosome.synthesize(rna_seq)
            if rna_seq:
                if protein_seq:  # if protein have something inside
                    seq_tuple = (srr_ret, rna_seq, protein_seq)
                else:
                    seq_tuple = (srr_ret, rna_seq, 'Non-coding RNA')
            else:
                seq_tuple = (srr_ret, 'Non-coding RNA')
            repertoire_list.append(seq_tuple)

        # return list of tuples for each genome
        return repertoire_list


class ProkaryoticCell(Cell):
    def __init__(self, genome, error_rate=0.0):
        self.name = 'ProkaryoticCell'
        self.num_copies = 1
        self.division_rate = 4
        self.genetic_code = {
            'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
            'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
            'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
            'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
            'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
            'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
            'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
            'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
            'UAC': 'Y', 'UAU': 'Y', 'UAA': None, 'UAG': None,
            'UGC': 'C', 'UGU': 'C', 'UGA': 'U', 'UGG': 'W'}
        self.start_codons = ['AUG', 'GUG', 'UUG']
        # __init__(self, name, genome, num_copies, genetic_code, start_codons, division_rate, error_rate=0.0)
        super().__init__(self.name, genome, self.num_copies, self.genetic_code, self.start_codons,
                         self.division_rate, error_rate)


class EukaryoticCell(Cell):
    def __init__(self, name, genome, division_rate, error_rate=0.0):
        self.name = name
        self.num_copies = 2
        self.start_codons = ['AUG']
        self.genetic_code = {'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
                             'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
                             'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
                             'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
                             'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
                             'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
                             'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
                             'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
                             'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
                             'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
                             'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
                             'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
                             'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
                             'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
                             'UAC': 'Y', 'UAU': 'Y', 'UAA': None, 'UAG': None,
                             'UGC': 'C', 'UGU': 'C', 'UGA': None, 'UGG': 'W'}
        self.division_rate = division_rate
        # __init__(self, name, genome, num_copies, genetic_code, start_codons, division_rate, error_rate=0.0)
        super().__init__(self.name, genome, self.num_copies, self.genetic_code, self.start_codons,
                         self.division_rate, error_rate)


class StemCell(EukaryoticCell):
    def __init__(self, genome, error_rate=0.0):
        self.division_rate = 3
        self.name = 'StemCell'
        # def __init__(self, name, genome, division_rate, error_rate=0.0)
        super().__init__(self.name, genome, self.division_rate, error_rate)


class NeuronCell(EukaryoticCell):
    def __init__(self, genome, error_rate=0.0):
        self.division_rate = 2
        self.name = 'NeuronCell'
        super().__init__(self.name, genome, self.division_rate, error_rate)


class MutantCell(StemCell):
    def __init__(self, genome, num_mutations):
        self.error_rate = 0.05
        super().__init__(genome, self.error_rate)
        self.name = 'MutantCell'
        self.num_mutations = num_mutations
        self.changed_seq = 0

    def mitosis(self):
        mitosis_list = []
        mutant_cell_divided = MutantCell(self.genome, self.num_mutations)
        new_genome = []
        for gene in mutant_cell_divided.genome:
            mutated_gene = mutant_cell_divided.cell_DNA_polymerase.transcribe(gene)
            mutant_cell_divided.num_mutations += mutant_cell_divided.cell_DNA_polymerase.num_of_errors
            # we add every time the number of errors that occurred
            mutant_cell_divided.cell_DNA_polymerase.error_rate = 0
            mutated_gene = mutant_cell_divided.cell_DNA_polymerase.transcribe(mutated_gene)
            # in order to get the original direction of "gene", with errors in it
            mutant_cell_divided.cell_DNA_polymerase.error_rate = 0.05
            # don't forget to return the number of errors to it's original rate
            new_genome.append(mutated_gene)
        mutant_cell_divided.genome = new_genome
        num_of_new_proteins = 0
        new_cell = StemCell(new_genome)
        for gene in new_genome:  # we check for number of mutated existing proteins and write them
            rna_seq = new_cell.cell_RNA_polymerase.transcribe(gene)
            protein_seq = new_cell.cell_ribosome.synthesize(rna_seq)
            if protein_seq:
                num_of_new_proteins += 1
        if mutant_cell_divided.num_mutations > 10:  # we convert the cell to cancer if 10 mutations occured
            mutant_cell_divided = CancerCell(mutant_cell_divided.genome, mutant_cell_divided.num_mutations)
        mutant_cell_divided.changed_seq = num_of_new_proteins
        for i in range(self.division_rate - 1):
            mitosis_list.append(mutant_cell_divided)
        self.changed_seq = 0  # we change the field to 0 to locate futher in the code new proteins
        mitosis_list.append(self)
        return mitosis_list

    def __mul__(self, other):
        new_divided_cells = []
        new_mutant_cell = MutantCell(self.genome, self.num_mutations)
        new_genome = []
        for gene in new_mutant_cell.genome:
            mutated_gene = new_mutant_cell.cell_DNA_polymerase.transcribe(gene)
            new_mutant_cell.num_mutations += new_mutant_cell.cell_DNA_polymerase.num_of_errors
            # we add every time the number of errors that occurred
            new_mutant_cell.cell_DNA_polymerase.error_rate = 0
            mutated_gene = new_mutant_cell.cell_DNA_polymerase.transcribe(mutated_gene)
            # in order to get the original direction of "gene", with errors in it
            new_mutant_cell.cell_DNA_polymerase.error_rate = 0.05
            # don't forget to return the number of errors to it's original rate
            new_genome.append(mutated_gene)
        new_mutant_cell.genome = new_genome
        if new_mutant_cell.num_mutations > 10:
            new_mutant_cell = CancerCell(new_mutant_cell.genome, new_mutant_cell.num_mutations)
        for number_of_divided in range(other):  # other is the number of cells replicated wrong
            new_divided_cells.append(new_mutant_cell)
        return new_divided_cells


class CancerCell(MutantCell):
    def __init__(self, genome, num_mutations=0):
        """

        :param genome: for mutant
        :param num_mutations: 0 default
        """
        super().__init__(genome, num_mutations)
        self.name = "CancerCell"
        self.division_rate = 10

    def mitosis(self):
        mitosis_list = []
        mutant_cell_divided = CancerCell(self.genome, self.num_mutations)
        new_genome = []
        for gene in mutant_cell_divided.genome:
            mutated_gene = mutant_cell_divided.cell_DNA_polymerase.transcribe(gene)
            mutant_cell_divided.num_mutations += mutant_cell_divided.cell_DNA_polymerase.num_of_errors
            # we add every time the number of errors that occurred
            mutant_cell_divided.cell_DNA_polymerase.error_rate = 0
            mutated_gene = mutant_cell_divided.cell_DNA_polymerase.transcribe(mutated_gene)
            # in order to get the original direction of "gene", with errors in it
            mutant_cell_divided.cell_DNA_polymerase.error_rate = 0.05
            # don't forget to return the number of errors to it's original rate
            new_genome.append(mutated_gene)
        mutant_cell_divided.genome = new_genome
        num_of_new_proteins = 0
        new_cell = StemCell(new_genome)  # what we do here is we check for mutated sequences, and write them down
        for gene in new_genome:
            rna_seq = new_cell.cell_RNA_polymerase.transcribe(gene)
            protein_seq = new_cell.cell_ribosome.synthesize(rna_seq)
            if protein_seq:
                num_of_new_proteins += 1
        mutant_cell_divided.changed_seq = num_of_new_proteins
        for i in range(self.division_rate - 1):  # we duplicate the same mutated cells
            mitosis_list.append(mutant_cell_divided)
        self.changed_seq = 0  # we change the field for further use for locating different seqs
        mitosis_list.append(self)
        return mitosis_list

    def __mul__(self, other):
        new_divided_cells = []
        new_mutant_cell = CancerCell(self.genome, self.num_mutations)
        new_genome = []
        for gene in new_mutant_cell.genome:
            mutated_gene = new_mutant_cell.cell_DNA_polymerase.transcribe(gene)
            new_mutant_cell.num_mutations += new_mutant_cell.cell_DNA_polymerase.num_of_errors
            # we add every time the number of errors that occurred
            new_mutant_cell.cell_DNA_polymerase.error_rate = 0
            mutated_gene = new_mutant_cell.cell_DNA_polymerase.transcribe(mutated_gene)
            # in order to get the original direction of "gene", with errors in it
            new_mutant_cell.cell_DNA_polymerase.error_rate = 0.05
            # don't forget to return the number of errors to it's original rate
            new_genome.append(mutated_gene)
        new_mutant_cell.genome = new_genome
        for number_of_divided in range(other):  # other is the number of cells replicated wrong
            new_divided_cells.append(new_mutant_cell)
        return new_divided_cells


class SequenceClassifier:
    def __init__(self, pattern_file):
        """
        pattern_file: input file.
        """
        self.dict_pattern_to_domain = self.__patterns_to_domains(pattern_file)

    def __convert_repetitions(self, protein, offset):
        """
        Converts expressions of (m, n) and (n) to their regex expressions
        :param protein: the protein
        :param offset: param to calc
        :return:
        """
        prosite_rep = protein[offset:]  # the part of the protein prosite we check for '()' brackets
        prosite_reps_len = len(prosite_rep)
        repetitions = ''
        if prosite_reps_len > 0:  # if there are brackets at all. if not, don't return a string
            if prosite_reps_len >= 3:  # if we have less then 3 characters it's for sure illegal prosite pattern
                if prosite_rep[0] == '(' and prosite_rep[prosite_reps_len - 1] == ')':  # check for legal brackets
                    repetitions += '{'  # replace '(' with regex '{'
                    pattern = '^[0-9]+,[0-9]+$|^[0-9]+$'  # the pattern to check for legal expresion inside brackets
                    checked_string = prosite_rep[1:prosite_reps_len - 1]  # the expresion inside brackets
                    if re.findall(pattern, checked_string):
                        # if check was successful it will have some list, a.k.a True
                        for index in range(1, prosite_reps_len - 1):  # -1 excludes the closing bracket
                            repetitions += prosite_rep[index]  # add the legal expresion to conversion
                    else:
                        return ValueError(protein)
                    repetitions += '}'  # replace ')' with '}' of regex
                else:
                    raise ValueError(protein)
            else:
                raise ValueError(protein)
        return repetitions

    def __check_what_amino_acids(self, prefix_bracket, protein, index_after_bracket):
        # amino_acids_letters = ['G', 'A', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S', 'P', 'V', 'I', 'C', 'Y', 'H',
        #                        'R', 'N', 'T', 'D']
        suffix_brackets = {'{': '}', '[': ']'}  # create a complementary dict for fewer code lines
        comp_bracket_index = None
        didnt_find_comp_bracket = True
        for letter_index in range(index_after_bracket, len(protein)):  # running inside the text inside the brackets
            if suffix_brackets[prefix_bracket] == protein[letter_index]:
                # check if the is a complementary bracket at all
                comp_bracket_index = letter_index  # remember it's index
                didnt_find_comp_bracket = False  # make know that we indeed found a complementary bracket
                break
        if didnt_find_comp_bracket:  # if we didn't find any bracket complementary, we raise error
            raise ValueError(protein)
        regex_check = protein[index_after_bracket: comp_bracket_index]  # string of text in brackets without them
        number_of_regex_repetitions = comp_bracket_index - index_after_bracket
        pattern = ("[GALMFWKQESPVICYHRNTD]" + "{" + str(number_of_regex_repetitions) + "}")
        if re.findall(pattern, regex_check):  # checking for legal amino acids inside brackets via findall method
            return regex_check
        else:
            raise ValueError(protein)

    def __switch(self, case, protein, index_of_case, index_of_protein, pattern_length):
        # check if the case argument matches the cases (usually the first letter in the protein pattern.
        # if '<' comes first then we check what's the first letter afterwards recursively
        # (if '<' comes again that's illegal and we check it))
        potential_end_of_protein = False
        if index_of_protein == pattern_length - 1:  # we always check whether the '>' comes in the right place in the
            # protein pattern. it can come only at the end of it
            potential_end_of_protein = True
        last_character = str(protein[len(protein) - 1])
        if case == "x":
            if last_character == '>' and potential_end_of_protein:
                return '.' + self.__convert_repetitions(protein[index_of_case:len(protein) - 1], index_of_case + 1) \
                       + '$'
            # convert x to ., check for '(n)' brackets of sorts and make it come in the end of seq via $
            else:
                return '.' + self.__convert_repetitions(protein, index_of_case + 1)
            # convert x to ., check for '(n)' brackets of sorts
        elif case == "{":
            regex_pattern = '[^' + self.__check_what_amino_acids('{', protein, index_of_case + 1) + ']'
            # convert {} to [^]
            comp_bracket_index = 3
            for letter_index in range(index_of_case + 1, len(protein)):  # finding '}' index
                if '}' == protein[letter_index]:
                    comp_bracket_index = letter_index
            if last_character == '>' and potential_end_of_protein:
                regex_pattern += self.__convert_repetitions(protein[0:len(protein) - 1], comp_bracket_index + 1) + '$'
                # check for '(n)' brackets of sorts and make it come in the end of seq via $
            else:
                regex_pattern += self.__convert_repetitions(protein, comp_bracket_index + 1)
                # check for '(n)' brackets of sorts
            return regex_pattern
        elif case == "[":
            regex_pattern = '[' + self.__check_what_amino_acids('[', protein, index_of_case + 1) + ']'
            # convert [] to []
            comp_bracket_index = 3
            for letter_index in range(index_of_case + 1, len(protein)):  # finding ']' index
                if ']' == protein[letter_index]:
                    comp_bracket_index = letter_index
            if last_character == '>' and potential_end_of_protein:
                regex_pattern += str(self.__convert_repetitions(protein[0:len(protein) - 1], comp_bracket_index + 1)) \
                                 + '$'
                # check for '(n)' brackets of sorts and make it come in the end of seq via $
            else:
                regex_pattern += str(self.__convert_repetitions(protein, comp_bracket_index + 1))
                # check for '(n)' brackets of sorts and make it come in the end of seq via $
            return regex_pattern
        elif case == "<":
            if index_of_protein != 0 or index_of_case != 0:  # it's relevant only if it comes
                # in first protein, and in the
                # beginning of it, or else it's an error
                raise ValueError(protein)
            regex_pattern = '^' + self.__switch(protein[index_of_case + 1], protein, index_of_case + 1,
                                                index_of_protein, pattern_length)
            # a recursive call for switch shouldn't get a stack overflow due
            # to it being called just once for code beauty
            return regex_pattern
        else:
            pattern = '^[GALMFWKQESPVICYHRNTD]{1}$'  # we check the case of a single amino acid starting the protein,
            # and
            # exclude any other case with an exception
            checked_string = case
            if re.findall(pattern, checked_string):
                if last_character == '>' and potential_end_of_protein:
                    return case + self.__convert_repetitions(protein[index_of_case:len(protein) - 1],
                                                             index_of_case + 1) + '$'
                    # write the amino acid, check for '(n)' brackets of sorts and make it come in the end of seq via $
                else:
                    return case + self.__convert_repetitions(protein, index_of_case + 1)
                    # write the amino acid, check for '(n)' brackets of sorts
            else:
                raise ValueError(protein)

    def prosite_to_python(self, pattern_dict):
        """
        :param: pattern_dict: the given dictionary the function will convert to python.
        :return: the updated dictionary by python rules.
        """
        regex_dict = {}  # dictionary of regex patterns to domains
        if pattern_dict:
            for prosite in pattern_dict:  # iteration every prosite pattern
                split_list_of_patterns = re.split('-', prosite)  # spliting the pattern via the '-'
                # which doesn't exist in regex
                regex_pattern = ''
                for protein_index in range(len(split_list_of_patterns)):  # iteraiting and converting
                    # every amino acid prosite
                    protein = split_list_of_patterns[protein_index]
                    regex_pattern += self.__switch(protein[0], protein, 0, protein_index, len(split_list_of_patterns))
                    # sendin
                    # the
                    # protein to be converted to regx
                regex_pattern = re.compile(regex_pattern)  # compiling the pattern strings we got
                regex_dict[regex_pattern] = pattern_dict[
                    prosite]  # replacing the old prosite with new regex patterns via
                # new dictionary of regex_dict
        return regex_dict

    def __patterns_to_domains(self, pattern_file):
        """

        :param pattern_file: path we need to check
        :return: dict of pattern (kay) to domain (value)
        """
        dic = {} #open new dict
        if pattern_file: # if there is input
            assert os.path.exists(pattern_file), pattern_file # check if valid path.
            if os.path.isfile(pattern_file): # check if file
                with open(pattern_file) as csv_file: #open csv (assume it is csv)
                    csv_reader = csv.reader(csv_file) # reader for lines run
                    for line in csv_reader: # for each line
                        if line[0] == 'Pattern': #title of csv , pass
                            pass
                        else: # a line that not a title
                            dic[line[0]] = line[1]
                # After the dictionary is done:
                return self.prosite_to_python(dic)
            else: # pattern file is not a valid file
                AssertionError(pattern_file)
        else:
            raise AssertionError("Not exist pattern_file path")

        return dic # return the final dictionary patterns_to_domains

    def classify(self, seq_list, csv_file):
        """
        :param: seq_list: list of sequences.
        :param: csv_file: csv file who will include the sequences after classification..
        """
        regex_dict = self.dict_pattern_to_domain  # regex to domains dictionary
        # open the file in the write mode
        with open(csv_file, 'w') as f:
            # create the csv writer
            w = csv.writer(f)
            # write the first row - titles to the csv file
            w.writerow(['Sequence', 'Domains'])
            protein_list = seq_list  # we get all the different proteins
            for protein in protein_list:  # locate every domain for every protein
                domains_string = ''
                for pattern in regex_dict:
                    if pattern.match(protein):  # we run every pattern on our protein
                        if domains_string:
                            domains_string += ';' + regex_dict[pattern]  # we check for various matches
                        else:
                            domains_string += regex_dict[pattern]
                if domains_string == '':  # we didnt get a single match
                    domains_string += 'NA'
                # we write the protein to domains match in our csv file
                w.writerow([protein, domains_string])


jsonCheck = {1: 'PatternsToDomainsInput', 2: 'GenomicSequencesInput', 3: 'OutputFile', 4: 'MaxCycles', 5: 'MaxCells'}


def main():
    random.seed(1)
    assert len(sys.argv) >= 2, "No argument inserted" #check if there are args
    assert len(sys.argv) == 2, "Only one input needed" #check if there are more then one arg
    msg = f"Not an exists path for file.json: {sys.argv[1]}"
    assert os.path.exists(sys.argv[1]), msg # check if path is exists
    jsonPath = sys.argv[1] # save it
    msg = f"Json file dont exists:{jsonPath}"
    assert os.path.isfile(jsonPath), msg # and finally check if there is a file in the path
    #we assume that when we got a file in the Path - the file is Json.
    try:
        # In any case, we'll make sure he's able to open the Jason, if the Jason is OK
        with open(jsonPath, 'r') as f:
            distros_dict = json.load(f) #load it
    except:
        raise AssertionError("Problem with open Json file to read from")
    msg = f"Json file must have 5 params. input file have {len(distros_dict.keys())} params"
    assert len(distros_dict.keys()) == 5, msg #check if there is 5 args.
    # param 1:
    msg = f"must have {jsonCheck[1]} key"
    assert jsonCheck[1] in distros_dict.keys(), msg # check if kay of first arg inside json
    msg = f"key without value {jsonCheck[1]}"
    assert distros_dict[jsonCheck[1]] != {}, msg # check if there is a value in this key
    msg = f"path to file.csv don't found: {distros_dict[jsonCheck[1]]}"
    assert os.path.isfile(distros_dict[jsonCheck[1]]), msg #check if path is to a file or not
    PatternsToDomainsInput = distros_dict[jsonCheck[1]]
    # param 2:
    msg = f"must have {jsonCheck[2]} key"
    assert jsonCheck[2] in distros_dict.keys(), msg # check if there is a key of second arg
    msg = f"key without value {jsonCheck[2]}"
    assert distros_dict[jsonCheck[2]] != {}, msg # check value
    msg = f"path to file.txt don't found: {distros_dict[jsonCheck[2]]}"
    assert os.path.isfile(distros_dict[jsonCheck[2]]), msg #check if file is valid
    GenomicSequencesInput = distros_dict[jsonCheck[2]] #save it
    # param 3:
    msg = f"must have {jsonCheck[3]} key"
    assert jsonCheck[3] in distros_dict.keys(), msg #check key
    msg = f"key without value {jsonCheck[3]}"
    assert distros_dict[jsonCheck[3]] != {} , msg # check value
    msg = f"path to file.csv don't found: {distros_dict[jsonCheck[3]]}"
    assert os.path.isfile(distros_dict[jsonCheck[3]]), msg #check id path is valid.
    OutputFile = distros_dict[jsonCheck[3]] # save it
    # param 4:
    # number of max divisions cycles - **include**
    msg = f"must have {jsonCheck[4]} key"
    assert jsonCheck[4] in distros_dict.keys(), msg #check key
    msg = f"key without value {jsonCheck[4]}"
    assert distros_dict[jsonCheck[4]] != {}, msg #check value
    MaxCycle = distros_dict[jsonCheck[4]]
    msg = f"value of {jsonCheck[4]} key is not an int"
    assert isinstance(MaxCycle, int), msg #check if num
    MaxCycle = int(MaxCycle) #save it
    assert MaxCycle >= 1, "MaxCycle range is 1-inf only."
    # param 5: number of max possible Cells
    msg = f"must have {jsonCheck[5]} key"
    assert jsonCheck[5] in distros_dict.keys(), msg #check key
    msg = f"key without value {jsonCheck[5]}"
    assert distros_dict[jsonCheck[5]] != {}, msg #check value
    MaxCells = distros_dict[jsonCheck[5]]
    msg = f"value of {jsonCheck[5]} key is not a number: {MaxCycle}"
    assert isinstance(MaxCells, int), msg #check if is int
    MaxCells = int(MaxCells) #Save it
    assert MaxCells > 1, "MaxCells range is 2-inf only."
    # print(PatternsToDomainsInput, GenomicSequencesInput, OutputFile, MaxCycle, MaxCells)
    try:
        #upolat the genome from txt file to a list []
        genomeMutantCell = [] #ths list for genome input
        with open(GenomicSequencesInput, 'r') as reader: #open read
            lines = reader.readlines() # lines the reader on them
            for line in lines:
                genomeMutantCell.append(line[:-1])
    except:
        raise AssertionError("Problem with uploading genome")
    try:
        #open new cell - MutantCell
        cell = MutantCell(genomeMutantCell, 0)
    except:
        raise AssertionError("Problem with initiate a Mutant Cell.")
    current_cells = [cell]  # the culture of cells in current
    limit_reached = False #boolien flag for not reaching limit
    sequence_classifier = SequenceClassifier(PatternsToDomainsInput)
    number_of_diff_proteins = 0  # counting the different number of proteins
    seq_list = []
    new_check_cell = StemCell(cell.genome)
    for gene in cell.genome:  # we check for coding seqs and add them to seq_list
        rna_seq = new_check_cell.cell_RNA_polymerase.transcribe(gene)
        protein_seq = new_check_cell.cell_ribosome.synthesize(rna_seq)
        if protein_seq:
            it_doesnt_contain_it = True
            for seq in seq_list:
                if seq == protein_seq:  # if the seq exists we dont add it
                    it_doesnt_contain_it = False
            if it_doesnt_contain_it:
                seq_list.append(protein_seq)
                number_of_diff_proteins += 1

    for i in range(int(MaxCycle)):  # we generate the maximal num of cycles
        start_num_cells = len(current_cells)
        for j in range(start_num_cells):  # run across this generation's cell
            next_cell = current_cells.pop(0)  # Remove the first cell from the list
            if len(current_cells) + next_cell.division_rate > int(MaxCells):  # if we can't division it
                # due to reaching max cells possible
                current_cells.append(next_cell)  # Just return the cell to the end of the list without replication
                limit_reached = True
                continue  # continue iterating, maybe we find a lesser division rate that can replicate
            next_gen = next_cell.mitosis()  # we can do mitosis
            current_cells += next_gen  # adding the next generation of mitosis to the cells
            for mutated_cell in next_gen:
                if mutated_cell.changed_seq != 0:  # checking for cells that mutated
                    number_of_diff_proteins += mutated_cell.changed_seq  # we count the changed proteins
                    new_check_cell = StemCell(cell.genome)  # create a cell for synthesizing
                    for gene in mutated_cell.genome:  # synthesize every gene in the mutated cell
                        rna_seq = new_check_cell.cell_RNA_polymerase.transcribe(gene)
                        protein_seq = new_check_cell.cell_ribosome.synthesize(rna_seq)
                        if protein_seq:  # if we synthesized something
                            it_doesnt_contain_it = True  # we hope it never appeared before
                            for seq in seq_list:
                                if seq == protein_seq:  # it appeared before, therefore we wont add it
                                    it_doesnt_contain_it = False
                            if it_doesnt_contain_it:
                                seq_list.append(protein_seq)  # add it to the seq list
                    break  # we do it once because all the mutations are equal
        if limit_reached and len(current_cells) == int(MaxCells):  # if we cant add division rate and reached
            # max cells possible, we dont run the program further
            break
    number_of_diff_proteins = len(seq_list)  # we fix the diff protein to be only the seqs that are in the list
    # because there are probably some repeating seqs
    print("Original cell: " + cell.__repr__())  # print first cell
    print("Final number of cells: " + str(len(current_cells)))  # print num of cells
    print('Protein repertoire size: ' + str(number_of_diff_proteins))  # print num of diffrent proteins
    max_num_of_mutations = 0
    for cell in current_cells:  # we extract the cells that mutated the most in current cells
        if cell.num_mutations > max_num_of_mutations:  # via this check
            max_num_of_mutations = cell.num_mutations  # update number of maximum mutations
    print('Mutations: ' + str(max_num_of_mutations))  # print num of mutations
    sequence_classifier.classify(seq_list, OutputFile)  # creating the CSV file of sequences to domains


if __name__ == "__main__":
    main()
