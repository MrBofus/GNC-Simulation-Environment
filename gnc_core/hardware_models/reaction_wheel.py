import numpy as np




def limit(val, maximum):
    if val > maximum:       return maximum
    elif val < -maximum:    return -maximum
    else:                   return val



class reactionWheelAssembly():
    def __init__(self, wheel_moi, max_wheel_speed, max_wheel_torque, min_wheel_torque):
        self.wheel_moi = wheel_moi
        self.max_wheel_speed = max_wheel_speed
        self.max_wheel_torque = max_wheel_torque
        
        self.wheel_speeds = np.array([0, 0, 0])
        self.wheel_torques = np.array([0, 0, 0])
        
        self.torque_command = [0, 0, 0]
    
    def commandReactionWheels(self, torque_command):
        t_command = []
        for i in range(3):
            t_command.append( limit(torque_command[i], self.max_wheel_torque) )
        
        self.torque_command = t_command
    
    def actuateReactionWheels(self, timestep):
        wheel_speeds_next = []
        for j in range(3):
            wheel_speed_candidate = self.wheel_speeds[j] + (timestep/self.wheel_moi)*self.torque_command[j]
        
            if abs(wheel_speed_candidate) > self.max_wheel_speed:
                wheel_speeds_next.append( self.wheel_speeds[j] )
            else:
                wheel_speeds_next.append( wheel_speed_candidate )
        
        wheel_torques_next = []
        for k in range(3):
            N = self.wheel_moi * (wheel_speeds_next[k] - self.wheel_speeds[k])/timestep
            wheel_torques_next.append( N )

        self.wheel_speeds = wheel_speeds_next
        self.wheel_torques = wheel_torques_next