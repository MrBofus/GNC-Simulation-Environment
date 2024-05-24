import numpy as np



def limit(val, maximum):
    if val > maximum:       return maximum
    elif val < -maximum:    return -maximum
    else:                   return val



class magnetorquerAssembly():
    def __init__(self, n_turns, length, width, max_current):
        self.n_turns = n_turns
        self.length = length
        self.width = width
        self.max_current = max_current
        
        self.current_command = np.array([0, 0, 0])
        self.torque = np.array([0, 0, 0])
    
    def commandMangetorquers(self, current_command):
        i_command = []
        for i in range(3):
            i_command.append( limit(current_command[i], self.max_current) )
        
        # i_command[2] = 0
        self.current_command = np.array(i_command)
    
    def actuateMagnetorquers(self, state):
        moment = self.n_turns * self.length * self.width * self.current_command
        self.torque = np.cross(moment, state.B_body)