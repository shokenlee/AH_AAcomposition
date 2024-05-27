#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
This modules contains the constants, functions and class for
AH search process based on the secondary structure results
'''

import math
import numpy as np


###################  Parameters  #######################

# Fixed parameters
amino_acids = ['A', 'G',  'V', 'L', 'I', 'F', 'W', 'M', 'Y', 'C',
       'S', 'T', 'R', 'K', 'N', 'Q', 'D', 'E', 'H', 'P']

hydrophobicity_dict = {'A': 0.31, 'R': -1.01, 'N': -0.6, 'D': -0.77, 'C': 1.54, 'Q': -0.22, 'E': -0.64, 'G': 0, 'H': 0.13,
                  'I': 1.8, 'L': 1.7, 'K': -0.99, 'M': 1.23, 'F': 1.79, 'P': 0.72, 'S': -0.04, 'T': 0.26, 'W': 2.25, 'Y': 0.96, 'V': 1.22}  # Based on Fauchere and Pliska 1983

charge_dict = {'A': 0, 'R': 1, 'N': 0, 'D': -1, 'C': 0, 'Q': 0, 'E': -1, 'G': 0, 'H': 0,
                  'I': 0, 'L': 0, 'K': 1, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0}

degree_step = 100  # Based on Eisenberg et al. 1982

minimum_dfactor = 0.68  # Based on HeliQuest help page


# Playable parameters
bulky_hydrophobic_residues = ['V', 'I', 'L', 'F', 'W', 'M', 'Y']
minium_helices = 6  # How many amino acids in both ends of an AA chunk have to be helix
minimum_len_consecutive_hydro = 3
minimum_hydrophobic_moment = 0.4


###################  Functions  #######################


def only_certain_characters(chunk, character_group):
    '''
    Judges if the given chunk contains a defined group of characters
    e.g. group = ['V', 'I', 'L', 'F', 'W']
         chunk = 'VLD' -> False; 'VLFWF' -> True
    e.g. group = 'H'
         chunk = 'HHHHHHH' -> True; 'HHHHHC' -> False
    '''
    
    for i in chunk:
        if i not in character_group:
            return False
    
    return True


def extract_protein_species_name(line):
    """
    Get the protein name and species name from the first line of the FASTA file
    """
    
    first_half = line[:line.find('OS=')]
    second_half = line[line.find('OS=') + 3:]
    protein_name = first_half[first_half.find(' '):].strip()
    species_name = second_half[:second_half.find('OX=')].strip()
    
    return protein_name, species_name



def get_face_seq(AA, vector_degree, degree_range, degree_step):
    """
    For a given sequence and given value of the degree (typically degree of the hydrophobic moment, 
    this function yields the sequence around the degree within the given range of degree (typically 90 or 45)
    when the sequence was plotted in a wheel projection
    
    e.g. "DAGMKRACGR", vector_degree = 175, degree_range = 90
    => face: "ARRGA", which is the sequence within 85 and 265 degrees in the wheel projection
    """
    
    face = ''
    
    # dictionary of AA index and its degree with a defined helix step (typically 100)
    # e.g. {0: 100, 1: 200, 2: 200, 3: 300, 4: 40, 5: 140, 6 : 240, 7: 340, 8: 80, ...}
    dic = {i: (degree_step * i) % 360 for i in range(len(AA))}
    
    # update the degree values in the dict around the vector such that
    # the vector is pointed toward 90 degree
    for k, v in dic.items():
        rotated_v = (v - vector_degree + 90) % 360
        dic[k] = rotated_v
    
    # sort the updated dictionary in order of the degree value
    # and retain those only with the degree range (e.g. from 0 to 180)
    keys = list(dic.keys())
    values = list(dic.values())
    dic_sorted = {keys[i]: values[i] for i in np.argsort(values) if (values[i] <= 90 + degree_range) & (values[i] >= 90 - degree_range)}
    
    # get the sequence of the face
    for k in dic_sorted.keys():
        face += AA[k]

    return face


###################  Class  #######################

        
class AA_seq():
    
    def __init__(self, sequence):
        self.sequence = sequence
        self.netcharge = None
        self.mean_hydrophobicity = None
        self.mean_hydrophobic_moment = None
        self.moment_degree = None
        self.dfactor = None
        self.hydro_phobic_face = None
        self.hydro_philic_face = None
        self.core_phobic_face = None
        self.contain_hydrophobic_cluster = None
        self.no_charge_hydrophobic_center = None
        
    def calculate_netcharge(self):
        
        # Prevent unitentional update 
        if self.netcharge == None:
            self.netcharge = 0
            for residue in self.sequence:
                self.netcharge += charge_dict[residue]
                
        return self.netcharge
    
    def calculate_hydrophobicity(self):
        
        # Prevent unitentional update
        if self.mean_hydrophobicity == None:
            total_hydrophobicity = 0
            for residue in self.sequence:
                total_hydrophobicity += hydrophobicity_dict[residue]
                self.mean_hydrophobicity = total_hydrophobicity / len(self.sequence)
        else:
            print("Seems like unintentional update is attempted")
                
        return self.mean_hydrophobicity
    

    def calculate_hydrophobic_moment(self):
        '''
        Return the magnitude of hydrophobic moment vector <µH> (mean_hydrophobic_moment)
        and its orientation
        e.g. 'FLNNAMSSLYSGWSSFTT' -> <µH> = 0.47.., degree = 90.1
        '''
        
        # Prevent unitentional update 
        if self.mean_hydrophobic_moment != None:
            print("Seems like unintentional update is attempted A")
        else:
            self.mean_hydrophobic_moment = 0

            x_hydrophobic_moment = 0
            y_hydrophobic_moment = 0
            
            # Calculate x, y of the hydrophobic moment vector
            for i, residue in enumerate(self.sequence):
                hydrophobicity = hydrophobicity_dict[residue]
                theta = degree_step * i
                x_hydrophobic_moment += hydrophobicity * math.cos(math.radians(theta))
                y_hydrophobic_moment += hydrophobicity * math.sin(math.radians(theta))

            hydrophobic_moment_len = math.sqrt(x_hydrophobic_moment **2 + y_hydrophobic_moment **2)  
            self.mean_hydrophobic_moment =  hydrophobic_moment_len / len(self.sequence) # Devided by the length by definition Eisenberg et al.
            
            # Calculate the degree
            if y_hydrophobic_moment >=0:
                self.moment_degree = math.acos(x_hydrophobic_moment / hydrophobic_moment_len) / math.pi * 180
            else:
                self.moment_degree = 360 - math.acos(x_hydrophobic_moment / hydrophobic_moment_len) / math.pi * 180
                
    
        
    def calculate_dfactor(self):
        '''
        Calculates d-factor based on mean hydrophobic moment and netcharge
        Prerequisite: mean_hydrohobic_moment has been executed
        '''
        # Prevent unitentional update 
        if self.dfactor != None:
            print("Seems like unintentional update is attempted B")
        else:
            self.dfactor = 0.944 * self.mean_hydrophobic_moment + 0.33 * self.calculate_netcharge()
            
        return self.dfactor

            
    def extract_face_sequences(self):
        '''
        Returns sequenes of hydrophobic and hydrophilic faces
        based on the given degree of the hydrophobic moment
        e.g. Input: 'RLRSAVSRAGSLLWMVAT' -> Output: 'AVALLVAGR' and 'WSTSSMRRL'
        
        It also gives a "core" or central hydrophobic face sequence.
        In the example above, ''
        '''
        
        # Prevent unintentional update
        if self.hydro_philic_face != None:
            print("Seems like unintentional update is attempted D")
        
        else:
            # get sequences of hydrophobic and hydrophilic faces and also "core" hydrophilic face
            # hydrophobic face is +/- 90 degrees from the hydrophobic moment
            # hydrophilic face is +/- 90 degrees from (the hydrophobic moment + 180 degrees) (ie the opposite)
            # core of the hydrophobic face is +/- 45 degrees from the hydrophobic moment
            
            self.hydro_phobic_face = get_face_seq(self.sequence, self.moment_degree, degree_range = 90, degree_step=degree_step)
            self.hydro_philic_face = get_face_seq(self.sequence, (self.moment_degree + 180) % 360, degree_range = 90, degree_step=degree_step)
            self.core_phobic_face = get_face_seq(self.sequence, self.moment_degree, degree_range = 45, degree_step=degree_step)
                
    
    def judge_hydrophobic_cluster(self, cluster_length = minimum_len_consecutive_hydro):
        '''
        Judge the presence of a cluster of bulky hydrophobic residues defined elsewhere
        The desired length of the cluster is typically 3
        e.g. 'RLRSAVSRAGSLFWMVAT' 
              -> hydrophobic face = 'AVALFVAGR', which contains 'LFV' 
                  -> True
        '''
        
        # Prerequisite: hydrophobic sequence has been obtained
        self.extract_face_sequences()
        
        # Prevent unintentional update
        if self.contain_hydrophobic_cluster != None:
            print("Seems like unintentional update is attempted E")
            
        else:
            self.contain_hydrophobic_cluster = False
            for i in range(len(self.hydro_phobic_face) - 2):
                if only_certain_characters(self.hydro_phobic_face[i:i+cluster_length], bulky_hydrophobic_residues):
                    self.contain_hydrophobic_cluster = True
                    break
                    
        return self.contain_hydrophobic_cluster
    
    def judge_no_charge_in_phobic_center(self):
        '''
        Judge the presence of a charged residue in the center of the hydrophobic face
        e.g. 'DAGLKRALGRRKGVWLRL' 
              -> center hydrophobic face = 'LRLWL', which contains 'R' 
                  -> False
        '''
        
        # Prerequisite: core hydrophobic sequence has been obtained
        if self.core_phobic_face == None:
            self.extract_face_sequences()
        
        # Prevent unintentional update
        if self.no_charge_hydrophobic_center != None:
            print("Seems like unintentional update is attempted F")
            
        else:
            self.no_charge_hydrophobic_center = True

            for charged_residue in [i for i in charge_dict.keys() if charge_dict[i] != 0]:
                if charged_residue in self.core_phobic_face:
                    self.no_charge_hydrophobic_center = False
                    break
        
        return self.no_charge_hydrophobic_center
        
   
    def judge_AH(self):
        """
        A given sequence is judged as a possible AH if it meets
        1) D factor > 0.68 and 2) hydrophobic cluster is present,
        which are extracted by the functions below
        """
        
        #
        self.calculate_hydrophobic_moment()
        
        # Consider it to be an AH candidate if
        # 1) D factor > 0.68 OR (µH > certain value (typically 0.4) & netcharge == 0)
        # 2) Hydrophobic cluster is present AND the center of hydrophobic face does NOT have charged residues
        if ((self.calculate_dfactor() > minimum_dfactor) | ((self.mean_hydrophobic_moment > minimum_hydrophobic_moment) & (self.netcharge == 0))) & (self.judge_hydrophobic_cluster() & self.judge_no_charge_in_phobic_center()):
            return True
        else:
            return False