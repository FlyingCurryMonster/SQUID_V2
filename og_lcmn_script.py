import logging
from matplotlib import units
import numpy as np

from sqlalchemy import Integer
from sympy import Float
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

import sys
sys.modules['cloudpickle'] = None  #THIS IS IMPORTANT TO INCLUDE IN CODE

from math import pi, sin, cos
import cmath
import time
from time import sleep
import pymeasure

from pymeasure.log import console_log

from pymeasure.instruments.srs import SR830
from pymeasure.instruments.tektronix import AFG3152C

from pymeasure.display.Qt import QtGui  
from pymeasure.display.windows import ManagedWindow
from pymeasure.display.widgets import PlotWidget, LogWidget, PlotWidget_datetimeaxis

from pymeasure.experiment import Procedure, Results, unique_filename
from pymeasure.experiment import IntegerParameter, FloatParameter, Parameter

from pyqtgraph import DateAxisItem

import pyvisa

def cubic(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d

unconstrain_poly_coeff = [-0.01440121, -0.00338684,  0.17276678,  0.04238221]
constrained_poly_coeff = [7.35395434E-6, -2.90946673E-2, 1.85983869E-1, 4.04657918E-2]

def chi_to_T(x, coeff = constrained_poly_coeff):
    Tsamp = np.arange(20, 0.33, step=-0.5E-3)
    inv_Tsamp = Tsamp**-1
    chi_fit = cubic(inv_Tsamp, *coeff)
    return np.interp(x, xp = chi_fit, fp = inv_Tsamp)**-1

class SQUID_Balancer_v2(Procedure): 
    
    delay = FloatParameter('delay', units = 's', default = 1)
    drive_V = FloatParameter('drive voltage', units = 'V', default = 3.25)
    cancel_gain = FloatParameter('cancellation gain', default = 1)
    squid_range = FloatParameter('SQUID output/Lockin input limit',units= 'V', default = 1.5E-3) #the max range on the lockin
    increment_voltage = FloatParameter('Increment voltage', units = 'V', default = 1E-3)
    rebalance_delay = IntegerParameter('Re-balance buffer count', default = 30)
    increment_sleep = FloatParameter('Increment sleep', units = 's', default = 5)

    sr830bus = IntegerParameter('SR830 bus', default = 2)
    sr830addr = IntegerParameter('SR830 address', default = 9)
    afgbus = IntegerParameter('AFG bus', default = 0)
    afgaddr = IntegerParameter('AFG address', default = 1)

    gain = FloatParameter('Tek/lock-in gain', default = 0.212)

    DATA_COLUMNS = ['UTC', 'timestamp', 'Temperature', 'tek_nx', 'tek_ny', 'X', 'Y', 'Tek amp', 'Tek phase', 'chi\'', 'chi\'\''] 
                
    sr830_full_address = 'GPIB' + str(sr830bus) + '::' + str(sr830addr) + '::' + 'INSTR' 
    afg_full_adress = 'GPIB' + str(afgbus) + '::' + str(afgaddr) + '::' + 'INSTR' 
    
    def find_gain(self):
        lockin = self.lockin_amp
        afg = self.afg
        gain = 1  #just some code to not cause an error 
        
        # zeros the input on the lock-in x channel with the tektronix and tells you what the gain is
        # gain = lockin.x/tek_amp
        self.amp = self.initial_amp
        afg.ch1.amp_vrms = self.amp
        afg.ch1.disable()  #set the amplitude to the lowest setting and disable to counterbalance
        signal = lockin.x  
        
        afg.ch1.enable()
        
        while abs(lockin.x) > 0.01*self.squid_range:
            self.tek_nx = self.tek_nx + 1
            self.phasor = complex(self.tek_nx*self.amp_increment, self.tek_ny*self.amp_increment)
            self.amp = self.initial_amp + abs(self.phasor)
            self.phase = self.initial_phase + cmath.phase(self.phasor)
            
            afg.ch1.amp_vrms = self.amp; #print(afg.ch1.amp_vrms)
            afg.ch1.phase_rad = self.phase

    def buffer_fill(self):
        for i in range(self.rebalance_delay + 1): 
            self.x_buffer.append(self.lockin_amp.x)
            self.y_buffer.append(self.lockin_amp.y)
            sleep(0.5)

    #log.info('SQUID_Balancer has been called')
    def startup(self):
        log.info('Starting SQUID_Balancer')
        self.x_buffer = []
        self.y_buffer = []

        self.lockin_amp = SR830(self.sr830_full_address)
        self.afg = AFG3152C(self.afg_full_adress)

        self.lockin_amp.sine_voltage = self.drive_V
        log.info('Drive voltage has been set to = {}'.format(self.lockin_amp.sine_voltage))
        self.chi_norm_constant = 5/self.lockin_amp.sine_voltage

        #self.x_buffer.append(self.lockin_amp.x)
        #self.y_buffer.append(self.lockin_amp.y)
        log.info('filling the buffer')
        self.buffer_fill()
        log.info('buffer filled')
        log.info('X buffer:')
        log.info(self.x_buffer)
        log.info('Y buffer:')
        log.info(self.y_buffer)

        self.tek_nx = 0
        self.tek_ny = 0
        self.amp_increment = self.increment_voltage
        #self.initial_amp = self.initial_tek_amp #7.1E-3
        self.initial_amp = self.afg.ch1.amp_vrms
        self.initial_phase = self.afg.ch1.phase_rad#.value * pi/180

        #self.phasor = complex(self.tek_nx*self.amp_increment, self.tek_ny*self.amp_increment)
        self.x_comp = self.initial_amp*cos(self.initial_phase) + self.tek_nx*self.amp_increment
        self.y_comp = self.initial_amp*sin(self.initial_phase) + self.tek_ny*self.amp_increment
        self.phasor = complex(self.x_comp, self.y_comp)

        self.amp = abs(self.phasor)
        self.phase = cmath.phase(self.phasor)

        self.afg.ch1.amp_vrms = self.amp
        self.afg.ch1.phase_rad = self.phase
        log.info('The initial amp has been set to %e' %self.afg.ch1.amp_vrms)
        log.info('The initial phase has been set to %e' %self.afg.ch1.phase_deg)

        self.t0 = time.time()
        
    def execute(self):
        lockin = self.lockin_amp
        afg = self.afg
        gain = self.gain
        
        while True:
            self.x_buffer.append(lockin.x)
            self.y_buffer.append(lockin.y)
            ###### ####### ####### ####### ####### ####### X CHANNEL ####### ####### ####### ####### ####### ####### ####### #######
            if all(abs(reading)>1.2*self.squid_range for reading in self.x_buffer[-1*self.rebalance_delay:]): #x_buffer used for x balancing
                log.info('Lockin X is out of squid_range')
                log.info(self.x_buffer[-1*self.rebalance_delay:])
                if lockin.x > 0:
                    log.info('Lockin X channel is large and positive')
                    while lockin.x > -0.8*self.squid_range:
                        self.x_buffer.append(lockin.x)
                        self.y_buffer.append(lockin.y)

                        self.tek_nx = self.tek_nx + 1
                        self.x_comp = self.initial_amp*cos(self.initial_phase) + self.tek_nx*self.amp_increment
                        self.y_comp = self.initial_amp*sin(self.initial_phase) + self.tek_ny*self.amp_increment
                        self.phasor = complex(self.x_comp, self.y_comp)
                        self.amp = abs(self.phasor)
                        self.phase = cmath.phase(self.phasor)
  
                        afg.ch1.amp_vrms = self.amp
                        afg.ch1.phase_rad = self.phase
                        
                        self.chi_re = (lockin.x + self.cancel_gain*gain*self.amp*cos(self.phase))*self.chi_norm_constant
                        self.chi_im = (lockin.y + self.cancel_gain*gain*self.amp*sin(self.phase))*self.chi_norm_constant

                        self.temp = chi_to_T(self.chi_re)
                        
                        data = {
                            'UTC': time.time() +2082844800, #the extra constant is to match Laview's timestamp start of 01/01/1904,
                            'timestamp': time.time()-self.t0,
                            'Temperature': self.temp,
                            'tek_nx' : self.tek_nx,
                            'tek_ny': self.tek_ny,
                            'X': lockin.x,
                            'Y': lockin.y,
                            'Tek amp': afg.ch1.amp_vrms,#self.amp,
                            'Tek phase': afg.ch1.phase_deg,#self.phase*180/pi,
                            'chi\'': self.chi_re,
                            'chi\'\'': self.chi_im 
                        }
                        self.emit('results', data)
                        log.debug('Lowering SQUID x output')
                        sleep(self.increment_sleep)

                        sleep(self.delay)
                          
                elif lockin.x < 0:
                    log.info('Lockin X channel is large and negative')
                    while lockin.x < 0.8*self.squid_range:
                        self.x_buffer.append(lockin.x)
                        self.y_buffer.append(lockin.y)

                        self.tek_nx = self.tek_nx-1
                        self.x_comp = self.initial_amp*cos(self.initial_phase) + self.tek_nx*self.amp_increment
                        self.y_comp = self.initial_amp*sin(self.initial_phase) + self.tek_ny*self.amp_increment
                        self.phasor = complex(self.x_comp, self.y_comp)
                        self.amp = abs(self.phasor)
                        self.phase = cmath.phase(self.phasor)
                        
                        afg.ch1.amp_vrms = self.amp
                        afg.ch1.phase_rad = self.phase
                        
                        self.chi_re = (lockin.x + self.cancel_gain*gain*self.amp*cos(self.phase))*self.chi_norm_constant
                        self.chi_im = (lockin.y + self.cancel_gain*gain*self.amp*sin(self.phase))*self.chi_norm_constant

                        self.temp = chi_to_T(self.chi_re)

                        data = {
                            'UTC' : time.time()+2082844800, #the extra constant is to match Laview's timestamp start of 01/01/1904,,
                            'timestamp': time.time()-self.t0,
                            'Temperature': self.temp,
                            'tek_nx' : self.tek_nx,
                            'tek_ny': self.tek_ny,
                            'X': lockin.x,
                            'Y': lockin.y,
                            'Tek amp': afg.ch1.amp_vrms,
                            'Tek phase': afg.ch1.phase_deg,
                            'chi\'': self.chi_re,
                            'chi\'\'': self.chi_im 
                        }
                        self.emit('results', data)
                        log.debug('Raising SQUID output')
                        sleep(self.increment_sleep)
                        sleep(self.delay)
                        
            ###### ####### ####### ####### ####### ####### Y CHANNEL ####### ####### ####### ####### ####### ####### ####### #######
            elif all(abs(reading) > 1.2*self.squid_range for reading in self.y_buffer[-1*self.rebalance_delay:]):
                 log.info('Lockin Y channel is out of squid_range')
                 log.info(self.y_buffer[-1*self.rebalance_delay:])
                 if lockin.y > 0:
                    log.info('Lockin Y channel is large and positive')
                    while lockin.x > -0.8*self.squid_range:
                        self.x_buffer.append(lockin.x)
                        self.y_buffer.append(lockin.y)

                        self.tek_ny = self.tek_ny + 1
                        self.x_comp = self.initial_amp*cos(self.initial_phase) + self.tek_nx*self.amp_increment
                        self.y_comp = self.initial_amp*sin(self.initial_phase) + self.tek_ny*self.amp_increment
                        self.phasor = complex(self.x_comp, self.y_comp)
                        self.amp = abs(self.phasor)
                        self.phase = cmath.phase(self.phasor)

                        
                        afg.ch1.amp_vrms = self.amp
                        afg.ch1.phase_rad = self.phase
                        
                        self.chi_re = (lockin.x + self.cancel_gain*gain*self.amp*cos(self.phase))*self.chi_norm_constant
                        self.chi_im = (lockin.y + self.cancel_gain*gain*self.amp*sin(self.phase))*self.chi_norm_constant

                        self.temp = chi_to_T(self.chi_re)

                        data = {
                            'UTC': time.time()+2082844800, #the extra constant is to match Laview's timestamp start of 01/01/1904,
                            'timestamp': time.time()-self.t0,
                            'Temperature': self.temp,
                            'tek_nx' : self.tek_nx,
                            'tek_ny': self.tek_ny,
                            'X': lockin.x,
                            'Y': lockin.y,
                            'Tek amp': afg.ch1.amp_vrms,
                            'Tek phase': afg.ch1.phase_deg,
                            'chi\'': self.chi_re,
                            'chi\'\'': self.chi_im 
                        }
                        self.emit('results', data)
                        log.debug('Lowering SQUID y output')
                        sleep(self.increment_sleep)

                        sleep(self.delay)
                          
                 elif lockin.y < 0:
                    log.info('Lockin Y channel is large and negative')
                    while lockin.x < 0.8*self.squid_range:
                        self.x_buffer.append(lockin.x)
                        self.y_buffer.append(lockin.y)

                        self.tek_ny = self.tek_ny-1
                        self.x_comp = self.initial_amp*cos(self.initial_phase) + self.tek_nx*self.amp_increment
                        self.y_comp = self.initial_amp*sin(self.initial_phase) + self.tek_ny*self.amp_increment
                        self.phasor = complex(self.x_comp, self.y_comp)
                        self.amp = abs(self.phasor)
                        self.phase = cmath.phase(self.phasor)
                        
                        afg.ch1.amp_vrms = self.amp
                        afg.ch1.phase_rad = self.phase
                        
                        self.chi_re = (lockin.x + self.cancel_gain*gain*self.amp*cos(self.phase))*self.chi_norm_constant
                        self.chi_im = (lockin.y + self.cancel_gain*gain*self.amp*sin(self.phase))*self.chi_norm_constant

                        self.temp = chi_to_T(self.chi_re)

                        data = {
                            'UTC' : time.time()+2082844800, #the extra constant is to match Laview's timestamp start of 01/01/1904,
                            'timestamp': time.time()-self.t0,
                            'Temperature': self.temp,
                            'tek_nx' : self.tek_nx,
                            'tek_ny': self.tek_ny,
                            'X': lockin.x,
                            'Y': lockin.y,
                            'Tek amp': afg.ch1.amp_vrms,
                            'Tek phase': afg.ch1.phase_deg,
                            'chi\'': self.chi_re,
                            'chi\'\'': self.chi_im 
                        }
                        self.emit('results', data)
                        log.debug('Raising y SQUID output')
                        sleep(self.increment_sleep)
                        sleep(self.delay)
                                       
            else:
                self.x_buffer.append(lockin.x)
                self.y_buffer.append(lockin.y)

                self.chi_re = (lockin.x + self.cancel_gain*gain*self.amp*cos(self.phase))*self.chi_norm_constant
                self.chi_im = (lockin.y + self.cancel_gain*gain*self.amp*sin(self.phase))*self.chi_norm_constant

                self.temp = chi_to_T(self.chi_re)

                data = {
                    'UTC': time.time()+2082844800, #the extra constant is to match Laview's timestamp start of 01/01/1904,
                    'timestamp': time.time()-self.t0,
                    'Temperature': self.temp,
                    'tek_nx' : self.tek_nx,
                    'tek_ny': self.tek_ny,
                    'X': lockin.x,
                    'Y': lockin.y,
                    'Tek amp': afg.ch1.amp_vrms,
                    'Tek phase': afg.ch1.phase_deg,
                    'chi\'': self.chi_re,
                    'chi\'\'': self.chi_im 
                }
                
                self.emit('results', data)
                log.debug('Lowering SQUID output')
                sleep(self.delay)
            
            if self.should_stop():
                log.warning('Caught the stop flag in the procedure')
                break
 
class LockinOutput(ManagedWindow): 
    def __init__(self):
        ChiPlot = PlotWidget_datetimeaxis(
            name = 'Chi', 
            columns = SQUID_Balancer_v2.DATA_COLUMNS,
            x_axis = 'UTC',
            y_axis = 'chi\'')

        #ChiPlot.plot.setAxisItems({'bottom': DateAxisItem(utcOffset= 2082844800)})

        TempPlot = PlotWidget_datetimeaxis(
            name = 'Temperature data', 
            columns = SQUID_Balancer_v2.DATA_COLUMNS, 
            x_axis = 'UTC', 
            y_axis = 'Temperature')

        inputs_displays = [ 
                     'delay',
                     'drive_V',
                     'cancel_gain',
                    'increment_voltage',
                     'squid_range',
                     'increment_sleep',
                     'rebalance_delay',
                     'gain',
                     'sr830bus',
                     'sr830addr',
                     'afgbus',
                     'afgaddr',
                     ]

        super().__init__(
            procedure_class=SQUID_Balancer_v2,
            inputs=  inputs_displays,
            displays = inputs_displays,
            x_axis = 'UTC',
            y_axis = 'X',
            directory_input= True,
            widget_list = (
                TempPlot, 
                ChiPlot,
            )
        )
        self.setWindowTitle('LCMN SQUID balancer')
        #D:\Data\LCMN SQUID
        self.directory = r'D:/Data/LCMN SQUID/'  
    
    def queue(self):
        directory = self.directory
        filename = unique_filename(directory, prefix = 'DATA')
        
        procedure = self.make_procedure()
        results = Results(procedure, filename)
        experiment = self.new_experiment(results)
        
        self.manager.queue(experiment)
    


if __name__ == '__main__':
    
    rm = pyvisa.ResourceManager()
    gpib_list = rm.list_resources()
    print(pymeasure.__version__)

    app = QtGui.QApplication(sys.argv)
    window = LockinOutput()
    window.show()
    sys.exit(app.exec_())