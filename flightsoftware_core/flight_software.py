import gnc_core.hardware_models.magnetometer as mm
import gnc_core.hardware_models.rate_gyro as rg
import gnc_core.hardware_models.star_tracker as st
import gnc_core.hardware_models.gps as gps

from flightsoftware_core.local_gnc_library import *




class schedulerApp():
    # measurements class {#4f7, 9}
    class measured_state():
        def __init__(self):
            self.measurements = {}

        def update(self, **kwargs):
            self.measurements = {}
            for arg, value in kwargs.items():
                self.measurements.update({arg:value})


    # init parameters {#121, 5}
    def __init__(self):
        self._p = {}
        self.m = self.measured_state()
        self.systemClock = 0


    # user-defined variables {#121, 31}
    def _update_user_variables(self, **kwargs):
        defaults = {
                    'bDot_gain':0.1,
            
                    'smc_kp':1.0,
                    'smc_kd':10.0,
                    'smc_sigma':0.1,
                    'smc_order':3,

                    'rg_noise':1/5000,
                    'rg_bias':0.001,

                    'st_noise':1/5000,
                    
                    'mg_noise':10**-11,
                    
                    'gps_noise_r':0.01,
                    'gps_noise_v':0.01,

                    'gs_range':800*10**3,

                    }
        
        self._p = defaults
        outstr = ''
        for arg, value in kwargs.items():
            if arg in defaults:
                outstr += '\t--- updated: ' + arg + '  ' + str(defaults[arg]) + ' --> ' + str(value) + '\n'
                self._p[arg]=value
            else:
                outstr += '\033[0;31m\t--- warning: ' + arg + ' not recognized as defualt variable, no action taken\033[0;32m\n'
        
        self.ground_station_list = [[30+41/60+27.8/3600, 88+10/60+31.5/3600, 0],
                                    [39.738868, -105.633692, 0],
                                    [4.811922, -55.840339, 0],
                                    [1.556200, 35.415574, 0],
                                    [50.454624, 8.355310, 0],
                                    [25.078205, 55.204769, 0],
                                    [23.303307, 89.085583, 0],
                                    [39.586561, 117.316221, 0],
                                    [1.310168, 103.700177, 0],
                                    [19.705257, -155.732732, 0],
                                    ]
        
        outstr += '\n' + str(self._p)
        printout(outstr)


    # interface functions {#121, 14}
    def _physics_to_harware_int(self, state):
        self._s = state
        self.systemClock += 1

    def _write_command_to_reactionWheel(self, command):
        self.rw_command = command
    
    def _write_command_to_magnetorquer(self, command):
        self.mt_command = command
    
    def _write_variables_to_state(self, mode, q_error):
        self.mode = mode
        self.q_error = q_error
    

    # main function {#112, 9}
    def _iterate(self):
        if self.systemClock == 1:
            self.mode = 'detumble'
        
        self._run_navigation()
        self._run_guidance()
        self._run_control()

        pass


    # flight apps {#211, 48}
    def _run_guidance(self):
        if self.mode == 'detumble':
            self.cmd = 'detumble'
            return None
        
        if qm.magnitude(checkVicinity(self.m.measurements, 
                                      self.ground_station_list, 
                                      self._p['gs_range'])) < self._p['gs_range']:
            self.cmd = 'downlink'
        
        else:
            self.cmd = 'nadir'


    def _run_navigation(self):
        wvec = rg.pull_gyro(self._s, noise=self._p['rg_noise'])
        qvec = st.pull_star_tracker(self._s, noise=self._p['st_noise'])
        bvec = mm.pull_magnetometer(self._s, noise=self._p['mg_noise'])

        rvec, vvec = gps.pull_gps(self._s, r_noise=self._p['gps_noise_r'], v_noise=self._p['gps_noise_v'])
        
        self.m.measurements.update(angularRate=wvec, quaternion=qvec, 
                                    bField=bvec, rvec=rvec, vvec=vvec)

    def _run_control(self):
        if self.mode == 'detumble':
            mt_command = bDotController(self.m.measurements['angularRate'], self.m.measurements['bField'], self._p['bDot_gain'])
            rw_command = [0, 0, 0]
            q_error = [0, 0, 0, 1]

            if qm.magnitude(self.m.measurements['angularRate']) < 0.002:
                self._write_command_to_magnetorquer(mt_command)
                self._write_command_to_reactionWheel(rw_command)
                self._write_variables_to_state('operation', q_error)
                return None
            
            self._write_command_to_magnetorquer(mt_command)
            self._write_command_to_reactionWheel(rw_command)
            self._write_variables_to_state(self.mode, q_error)
            return None
        
        elif self.mode == 'operation':
            if self.cmd == 'nadir':
                q_setpoint = ground_target_guidance(self.m.measurements, 'nadir', 
                                                    self.ground_station_list, self._p['gs_range'])
   
            elif self.cmd == 'prograde':
                q_setpoint = ground_target_guidance(self.m.measurements, 'prograde', 
                                                    self.ground_station_list, self._p['gs_range'])
            
            elif self.cmd == 'downlink':
                q_setpoint = ground_target_guidance(self.m.measurements, 'downlink', 
                                                    self.ground_station_list, self._p['gs_range'])
            
            else:
                q_setpoint = [-1, -1, -1, -1]

            
            q_error, rw_command = slidingModeController(self.m.measurements['angularRate'], 
                                                        self.m.measurements['quaternion'], q_setpoint, 
                                                        self._p['smc_kp'], self._p['smc_kd'], self._p['smc_sigma'], self._p['smc_order'])
            
            '''
            rw_command = [0, 0, 0]
            # mt_command = [0, 0, 0]
            m_gain = 3*10**3
            q_error, mt_command = magneticController( self.m.measurements['angularRate'], 
                                                      self.m.measurements['quaternion'], [0, 0, 0, 1], 
                                                      m_gain*0.05, m_gain*0.1, self._p['smc_sigma'], self._p['smc_order'],
                                                      self.m.measurements['bField'] )
            '''
            
            mt_command = [0, 0, 0]
            self._write_command_to_magnetorquer(mt_command)
            self._write_command_to_reactionWheel(rw_command)
            self._write_variables_to_state(self.mode, q_error)
            return None