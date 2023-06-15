import numpy as np


"""""
type of state functions :
is reachable      (Done)
is absorbing     (Done)
is closed set   (canceled)
is not irreducable  (done)
is recurrent     (Done)
is transient     (Done)
*************************************
type of markov chain class functions :
is reducible   (done)
is priodic/ aperiodic     (Done)
is absorbing (Absorbing if every state can reach an absorbing state.)     (done)
is ergodicity   (done)


"""
class MarkovChain:
    def __init__(self, transition_matrix, initial_probabilities):
        self.transition_matrix = transition_matrix
        self.initial_probabilities = initial_probabilities
        self.num_states = len(transition_matrix)
        self.absorbing_state_ID=None
        self.peridic,self.recurrent=False,False
        self.irreducible=False

    def is_absorbing_state(self,state_ID):
        # an absorbing state is a state where the probability of being in the state is 100%
        # that means if the corresponding element in diagonal in the transition matrix is 1 , the state is absorbing
        if self.transition_matrix[state_ID][state_ID]==1:
                self.absorbing_state_ID=state_ID
                return "is absorbing state"
        else:
            return "is not absorbing state"
    
    def is_not_reachable_state(self,state_ID):
        # the state is not reachable if all values in its column in the transient matrix are 0%
        if self.transition_matrix[:][state_ID].sum()==0:
                self.reducible=True
                return "is not reachable state"
        return "is reachable state"

    def is_transient_state(self, state_ID):
        # the state is transient if the probability to reach it is less then 100%
        if np.all(self.transition_matrix[state_ID, :] == 0):
            return "is transient state"
        else:
            self.recurrent=True
            return "is recurrent state"
    
    def Type_Of_State(self, state_id):
        # determine the type of the state with specific ID
        print("the state {} ".format(state_id),self.is_absorbing_state(state_id))
        print("the state {} ".format(state_id),self.is_not_reachable_state(state_id))
        print("the state {} ".format(state_id),self.is_transient_state(state_id))
    
    def is_periodic(self):
        # the type of markov chain class is periodic if the greatest common divisor d of the set is ≥ 2 (regular time), else it is aperiodic
        eigenvalues, _ = np.linalg.eig(self.transition_matrix)
        if not np.all(np.isclose(np.abs(eigenvalues), 1)) :
            self.peridic=True
            return "is periodic"
        else:
            return "is aperiodic"

    def Type_Of_Markov_Chain_Class(self):
        print("the markov chain class ",self.is_periodic())
        if self.is_not_reachable_state(self.absorbing_state_ID)=="is reachable state":
            print("the markov chain class is absorbing")
        if self.peridic and self.recurrent :
            print("the markov chain class is ergodicity")
        if self.irreducible:
            print("the markov chain class is irreducible")
        else:
            print("the markov chain class is not irreducible")

    def get_steady_state_probabilities(self):
        
        eigenvalues, eigenvectors = np.linalg.eig(self.transition_matrix.T)
        steady_state_index = np.where(np.isclose(eigenvalues, 1))[0]
        steady_state_vector = eigenvectors[:, steady_state_index].real.flatten()
        steady_state_vector /= np.sum(steady_state_vector)
        return steady_state_vector
        

    def get_state_probabilities_after_steps(self, num_of_steps):
        state_probabilities_after_steps= np.dot(self.initial_probabilities,self.transition_matrix)
        for i in range(num_of_steps-1):
          state_probabilities_after_steps = np.dot(state_probabilities_after_steps,self.transition_matrix)
        return state_probabilities_after_steps


# take inputs
num_of_states=int(input("Enter the number of states : "))
transition_matrix=np.zeros([num_of_states,num_of_states],float)
for i in range(num_of_states):
   for j in range(num_of_states):
      transition_matrix[i][j]=float(input("Enter the transition probability from state {} to state {} : ".format(i,j)))

initial_probabilities=np.zeros(num_of_states,float)
for i in range(num_of_states): 
  initial_probabilities[i]= float(input("Enter the initial probabilities for state {} : ".format(i)))
num_of_steps=int(input("Enter the number of steps to calculate the state probabilities after them : "))
markov_chain = MarkovChain(transition_matrix, initial_probabilities)

print("\n\n the results are : \n")
# type of each state 
for i in range(num_of_states):
 markov_chain.Type_Of_State(i)

# type of markov chain class
markov_chain.Type_Of_Markov_Chain_Class()

# calculate the steady state probabilities and the state probavilities after number of steps 
print("the steady state probabilities : ",markov_chain.get_steady_state_probabilities(),"\nthe state probabilities after {} steps".format(num_of_steps),markov_chain.get_state_probabilities_after_steps(num_of_steps))
