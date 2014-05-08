# (C) 2012-2013 Elke Schaper

import os, shutil, subprocess, tempfile
import numpy as np
import logging
logger = logging.getLogger('root')

from ..repeat import repeat_io
#from ..repeat.repeat_score import loadModel
from . import hmm_io

################################### HMM class #########################################

class HMM:
    """ A cyclic HMM applicable to describe sequence Tandem Repeats """
    
    def __init__(self, tandem_repeat=None, hmm_file = None, accession = None, hmm_file_copy = None, id = None):
        
        if id:
            self.id = id
        else:
            self.id = "unnamed"
        
        if hmm_file and accession:
            hmmer_probabilities = hmm_io.read_HMMER(hmm_file, id = accession)
            if not hmmer_probabilities:
                print("<hmmer_probabilities> were not found :(")
            else:
                self.HMM_from_model(hmmer_probabilities)
        elif tandem_repeat:
            self.HMM_from_TR(tandem_repeat, hmm_file_copy = hmm_file_copy)
        else:
            self.HMM_example()
        
    def HMM_example(self):    
        #states = ["N", "B", "M1", "M2", "M3", "E", "C"]
        self.states = ["N", "M1", "M2", "C"]
                
        ## Initialisation
        self.p_t = { iS: { iS2: 0  for iS2 in self.states}  for iS in self.states[:-1]}
        self.p_0 = { iS: 1/3  for iS in self.states }
        
        ## Transition probabilities
        # Feed Values to p_t
        self.p_t["N"] = {"N": 0.5, "M1": 0.5}
        self.p_t["M1"] = {"M2": 0.5, "C": 0.5}
        self.p_t["M2"] = {"M1": 0.5, "C": 0.5}
        self.p_t["C"] = {"C": 1}
        
        # emissions
        self.emissions = ["A", "C", "G", "T"]
        
        # emission probabilities
        self.p_e = { iS: { iE: 0.25  for iE in self.emissions  }  for iS in self.states }
        self.p_e['M1'] = {"A": 0.9, "C": 0.025, "G": 0.025, "T": 0.025}
        self.p_e['M2'] = {"A": 0.025, "C": 0.9, "G": 0.025, "T": 0.025}
     
    def HMM_from_TR(self, tandem_repeat, hmm_file_copy = None):
        """ Let HMMbuild build the HMM from <tandem_repeat>
            Next, create circular HMM_from_model(). """
    
        # Create a temporary directory
        tmp_dir = tempfile.mkdtemp()
        
        # Save TR as Stockholm file
        stockholm_file = os.path.join(tmp_dir, self.id+".sto")
        repeat_io.save_repeat_stockholm(tandem_repeat.msaD, stockholm_file)
        
        # Run HMMbuild to build a HMM model, and read model
        p = subprocess.Popen(["hmmbuild", "--amino", "tmp.hmm", self.id+".sto"], stdout=subprocess.PIPE, stderr=None, cwd=tmp_dir)
        p.wait()
        
        hmm_file = os.path.join(tmp_dir, "tmp.hmm")
        if hmm_file_copy:
            shutil.copy(hmm_file, hmm_file_copy)
        hmmer_probabilities = hmm_io.read_HMMER(hmm_file, id = self.id, type = "NAME")
        
        shutil.rmtree(tmp_dir)
                
        self.HMM_from_model(hmmer_probabilities)
              
    def HMM_from_model(self, hmmer_probabilities):
    
        """ Load a HMM from hmm_file
            Store probabilities as log10arithms.
            
            Parameters:
            <hmm_file> = 'path/to/file.hmm'
            <accession> = 'PF00560'
            <prior_indel_insertion> = {'mu': 0.5, 'sigma_squared': 0.81}
            
            
        """
                
        self.hmmer = hmmer_probabilities
        self.alphabet = self.hmmer['letters']
    
        # Warning! The following line is only appropriate, if there are only "letters", "COMPO", and the match state keys in <hmmer_probabilities>
        self.lD = len(hmmer_probabilities.keys()) - 2
        
        # Assume: sequence_type = 'AA'
        #Q,null_model_emission_p,alphabet = loadModel('lg')
  
        # Initialise all HMM states to default value (e.g. transition to or from terminal state all have the same cost 0).
        self.initialise_HMM_structure(self.lD)

        dTranslate_States = {i: i[1:] for i in self.match_states+self.insertion_states}
        for i in self.terminal_states:
            dTranslate_States[i] = "COMPO"
        
        # Store null model emission probabilities for ["N","C"] from hmmer_probabilities["COMPO"]
        # Store match_state and insertion_state emission probabilities for ["M0",...] from hmmer_probabilities["0"],...
        for iState in self.match_states:
            iEmission_Probabilities = self.hmmer[dTranslate_States[iState]]['emissions'][:len(self.alphabet)]
            self.set_emission_probability_hmmer3(iState, iEmission_Probabilities)
        for iState in (self.insertion_states + self.terminal_states):
            iEmission_Probabilities = self.hmmer[dTranslate_States[iState]]['insertion_emissions'][:len(self.alphabet)]
            self.set_emission_probability_hmmer3(iState, iEmission_Probabilities)
            
        # Store transition probabilities (m->m    m->i    m->d    i->m    i->i    d->m    d->d):
        # First: All state defined in HMMer model
        for i,iState in enumerate(self.match_states):
            iTransition_Probabilities = self.hmmer[dTranslate_States[iState]]['transition']
            # Completing the circle of the HMM.
            if i == self.lD - 1:
                iTransition_Probabilities = self.set_circle_transition_probability_hmmer3(iTransition_Probabilities, dTranslate_States[iState])
            self.set_transition_probability_hmmer3(i, iTransition_Probabilities)
            
        # Translate deletion states into direct transitions from match to match states.
        p_t_d = {}
        for i,iMatch_state in enumerate(self.match_states):
             p_t_d[iMatch_state] = self.get_direct_transition_probabilities_for_deletions(i,iMatch_state)
             #print(iMatch_state)
             #print('M' + str(((i+1)%self.lD)+1))
             #print(p_t_d[iMatch_state]['M' + str(((i+1)%self.lD)+1)])
        for iMatch_state in self.match_states:
            for iGoal_state, iP in p_t_d[iMatch_state].items():
                self.p_t[iMatch_state][iGoal_state] = iP
        
        # As it is never likely to got from a terminal state to a deletion state or vice versa, this transition probabilities can be ignored.
        # Thus, you may advance and:
        # Delete everything you ever knew about deletion states.
            
        self.states = [self.terminal_states[0]] + self.match_states + self.insertion_states + [self.terminal_states[1]]
        self.p_t = {iState: {i:j for i,j in self.p_t[iState].items() if not i in self.deletion_states} for iState in self.states}
        
        del self.deletion_states
        
        
    def initialise_HMM_structure(self,lD):
        
        '''
        Initialise all states
        Set transition probabilities to None/0
        Set initial probability for all states equal.
        '''
    
        ## Build a HMM from these likelihoods.
        self.match_states = ["M{0}".format(str(i+1)) for i in range(lD)]
        self.insertion_states = ["I{0}".format(str(i+1)) for i in range(lD)]
        self.deletion_states = ["D{0}".format(str(i+1)) for i in range(lD)]
        self.terminal_states = ["N","C"]
        self.states = [self.terminal_states[0]] + self.match_states + self.insertion_states + self.deletion_states + [self.terminal_states[1]]
        logger.debug("HMM states: {0}".format(self.states))      
                
        ## Initialisation
        # The transition probability is initially set to None for all states.
        self.p_t = { iS: { iS2: None  for iS2 in self.states}  for iS in self.states}
        # The initial probability is equal in all states. (np.log10(1)=0)
        self.p_0 = { iS: 0  for iS in self.states }
                
        ## Transition probabilities
        
        # First, mark which transitions are not penalised, e.g. bear a 'probability of np.log10(1)=0':
        # - Transitions from "N" to any "M_k" 
        # - Transitions from any "M_k" to "C"
        # - Transitions from "M_k" to "M_k+1" (for k<n) and "M_k" to "M_0" (for k==n)
        # - Transitions from "N" to "N"
        # - Transitions from "C" to "C"
        for i in range(lD):
            self.p_t["N"]["M{0}".format(i%lD + 1)] = 0
            self.p_t["M{0}".format(i%lD + 1)]["C"] = 0
            self.p_t["M{0}".format(i%lD + 1)]["M{0}".format((i+1)%lD + 1)] = 0
        self.p_t["N"]["N"] = 0
        self.p_t["C"]["C"] = 0
        
        self.p_e = {}
        
    def set_emission_probability_hmmer3(self, state, lEmission_Probabilities):
        '''
        Set p_e of states to <lEmission_Probabilities> for <state> given <self.alphabet>
        In HMMER3 data, emission probabilities are -ln(p).
        Return log10(p), i.d. convert between the two. Conversion: p_Local = - p_HMM * log10(e)
    
        Parameters (e.g.):
        state = 'M4'
        lEmission_Probabilities = ['3.27687', '2.31397', '3.32252', '3.12746', '2.89175', '3.34719', '2.28730', '3.54139', '2.53154', '2.64774', '3.75733', '3.23860', '3.57894', '3.26290', '2.66343', '2.61544', '2.91770', '3.26739', '5.09378', '3.16816']
        self.alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
        Return format (pseudo code):
        [{iA: np.log10(p(iA,iM)) for iA in alphabet.keys()} for iM in lMatch]

        '''
        lEmission_Probabilities = [-float(i)*np.log10(np.exp(1)) for i in lEmission_Probabilities]
        self.p_e[state] = {iL: iEP for iL,iEP in zip(self.alphabet, lEmission_Probabilities)}


    def set_circle_transition_probability_hmmer3(self, lTransition_Probabilities, final_state):
        ''' In HMMER3 models, no transition in circle from the last match state to the first
        match state are included.
        
        Values need to be found for the following transitions:
        m->m    m->i    m->d    i->m    i->i    d->m    d->d
        
        Compare ftp://selab.janelia.org/pub/software/hmmer3/3.0/Userguide.pdf pp.84-85.
        
        Input: -ln(p) Output: Also -ln(p')
        '''
        
        # M->M, M->I, M->D: Use an average from all other transitions from match states.
        mean_M_M = np.mean([np.exp(-float(iP["transition"][0])) for iState, iP in self.hmmer.items() if iState not in ['COMPO', 'letters', final_state]])
        mean_M_I = np.mean([np.exp(-float(iP["transition"][1])) for iState, iP in self.hmmer.items() if iState not in ['COMPO', 'letters', final_state]])
        mean_M_D = np.mean([np.exp(-float(iP["transition"][2])) for iState, iP in self.hmmer.items() if iState not in ['COMPO', 'letters', final_state]])
        sum_p = np.sum([mean_M_M,mean_M_I,mean_M_D])
        lTransition_Probabilities[:3] = [-np.log(i/sum_p) for i in [mean_M_M,mean_M_I,mean_M_D]]

        # I->M and I->I: Use an average from all other transitions from insertion states.
        mean_I_M = np.mean([np.exp(-float(iP["transition"][3])) for iState, iP in self.hmmer.items() if iState not in ['COMPO', 'letters', final_state]])
        mean_I_I = np.mean([np.exp(-float(iP["transition"][4])) for iState, iP in self.hmmer.items() if iState not in ['COMPO', 'letters', final_state]])
        sum_p = np.sum([mean_I_M,mean_I_I])
        lTransition_Probabilities[3:5] = [-np.log(i/sum_p) for i in [mean_I_M,mean_I_I]]
        
        # D->M and D->D: Use an average from all other transitions from deletion states.
        mean_D_M = np.mean([np.exp(-float(iP["transition"][5])) for iState, iP in self.hmmer.items() if iState not in ['COMPO', 'letters', final_state]])
        mean_D_D = np.mean([np.exp(-float(iP["transition"][6])) for iState, iP in self.hmmer.items() if iState not in ['COMPO', 'letters', final_state]])
        sum_p = np.sum([mean_D_M,mean_D_D])
        lTransition_Probabilities[5:] = [-np.log(i/sum_p) for i in [mean_D_M,mean_D_D]]
        
        return lTransition_Probabilities
        
    def set_transition_probability_hmmer3(self, state_index, lTransition_Probabilities):
        '''
        Set p_t of states to <lTransition_Probabilities> for the <state_index>th state
        In HMMER3 data, transition probabilities are -ln(p).
        Return log10(p), i.d. convert between the two. Conversion: p_Local = - p_HMM * log10(e)
    
        Parameters (e.g.):
        state = 4
        lTransition_Probabilities = ['0.00021', '8.85487', '9.57722', '0.61958', '0.77255', '0.48576', '0.95510']

        '''
        lTransition_Probabilities = [-float(i)*np.log10(np.exp(1)) for i in lTransition_Probabilities]
        
        # Match state state_index -> Match state state_index + 1
        iM = self.match_states[state_index]
        iM_1 = self.match_states[(state_index + 1)%self.lD]
        self.p_t[iM][iM_1] = lTransition_Probabilities[0]
        
        # Match state state_index -> Insertion state state_index
        self.p_t[iM][self.insertion_states[state_index]] = lTransition_Probabilities[1]
        
        # Match state state_index -> Deletion state state_index + 1
        iD_1 = self.deletion_states[(state_index + 1)%self.lD]
        self.p_t[iM][iD_1] = lTransition_Probabilities[2]
        
        # Insertion state state_index ->  Match state state_index + 1  
        iI = self.insertion_states[state_index]
        self.p_t[iI][iM_1] = lTransition_Probabilities[3]

        # Insertion state state_index ->  Insertion state state_index
        iI = self.insertion_states[state_index]
        self.p_t[iI][iI] = lTransition_Probabilities[4]
        
        # Deletion state state_index -> Match state state_index + 1
        iD = self.deletion_states[state_index]
        self.p_t[iD][iM_1] = lTransition_Probabilities[5]

        # Deletion state state_index -> Deletion state state_index + 1
        self.p_t[iD][iD_1] = lTransition_Probabilities[6]
                
    def get_direct_transition_probabilities_for_deletions(self,state_index,state):
        ''' Translate all deletions following <state> with <state_index>
            into direct match state to match state transitions.
            Probabilities need to be multiplied along the edges (log(p) need to be summed).
            
            Return transition probabilities
        '''
        
        transition = {}
        p_Deletion_State = self.p_t[state][self.deletion_states[(state_index + 1)%self.lD]]
        for iState_index in range(1,self.lD):
            iDeleted_state_index = (state_index + iState_index)%self.lD
            iDeleted_state = self.deletion_states[iDeleted_state_index]
            iGoal_state_index = (state_index + iState_index + 1)%self.lD
            iGoal_state = self.match_states[iGoal_state_index]
            transition[iGoal_state] = p_Deletion_State + self.p_t[iDeleted_state][iGoal_state]
            p_Deletion_State += self.p_t[iDeleted_state][self.deletion_states[iGoal_state_index]]
        
        return transition
        
    
############################### External parameters ######################################


def hmmer3_emission_probabilities(hmmer_probabilities, letters, lMatch):

    '''
    Get emission probabilities from hmmer3 hmm file.
    In hmm file, emission probabilities are -ln(p).
    Return log10(p), i.d. convert between the two. Conversion: p_Local = - p_HMM * log10(e)
    
    Parameters (e.g.):
    letters = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']
    lMatch = ['M'+str(i) for i in range(24)]
    
    Return format (pseudo code):
    [{iA: np.log10(p(iA,iM)) for iA in alphabet.keys()} for iM in lMatch]

    '''
    
    # Test: Is the number of match states in both models equal?
    if not len(hmmer_probabilities.keys())-2 == len(lMatch):
        print('Match states HMMER: {0} Match states local: {1}'.format(len(hmmer_probabilities.keys())-2, len(lMatch)))
        raise ValueError('The number of match states in HMMer model and local model does not match')
        
    # Test: Are all <letters> represented in the HMMER HMM?
    if any( iL not in hmmer_probabilities['letters'] for iL in letters ):
        missing = [iL for iL in letters if iL not in hmmer_probabilities['letters']]
        print('Missing representation in Hmmer File: {0}'.format(missing))
        raise ValueError('Some letters in the local HMM are not represented in the HMMER HMM.')
    
    return [{iL: -float(iP)*np.log10(np.exp(1)) for iL,iP in zip(hmmer_probabilities['letters'], data['emissions']) if iL in letters} for key,data in hmmer_probabilities.items() if key not in ['letters','COMPO']]
    
     
##################################### Tests ##############################################

     
def test():

    ''' To be implemented... '''    
    tandem_repeat = ...
    divergence = 0
    calculate_log10_offspring_likelihood(tandem_repeat, divergence) 
     
##################################### Main ###############################################

def main():
    my_TR = repeat_info.Repeat(begin = 0, msa = ['A-G', 'ACG', 'ACG'], sequence_type = 'DNA')
    #my_TR = TR()
    my_HMM = HMM(my_TR, divergence = 0.1)
    
    
if __name__ == "__main__":  
    main()