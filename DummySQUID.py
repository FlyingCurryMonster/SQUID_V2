import logging
import numpy as np
from utils import FixedSizeBuffer

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

import sys

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
    # drive_V = FloatParameter('drive voltage', units = 'V', default = 3.25)
    cancel_gain = FloatParameter('cancellation gain', default = 1)
    squid_range = FloatParameter('SQUID output/Lockin input limit',units= 'V', default = 1.5E-3) #the max range on the lockin
    increment_voltage = FloatParameter('Increment voltage', units = 'V', default = 1E-3)
    buffer_size = IntegerParameter('Re-balance buffer count', default = 30)
    increment_sleep = FloatParameter('Increment sleep', units = 's', default = 5)


    sr830addr = Parameter('SR830 address', default='2::9')
    afgaddr = Parameter('AFG address', default='0::1')

    teksr_gain = FloatParameter('Tek/lock-in gain', default = 0.212)

    inputs_displays = [ 
                    'delay',
                    'drive_V',
                    'cancel_gain',
                    'increment_voltage',
                    'squid_range',
                    'increment_sleep',
                    'buffer_size',
                    'teksr_gain',
                    'sr830addr',
                    'afgaddr',
                    ]

    DATA_COLUMNS = [
        'UTC', 'timestamp',
        'lockin_drive', 'X', 'Y',
        'tek_nx', 'tek_ny',
        'X_cancel', 'Y_cancel', 
        'chi_real', 'chi_imag'
        ] 
                
    def find_gain(self):
        lockin = self.lockin
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
    
    # USE FIXED SIZE BUFFER METHODOLOGY
    def buffer_fill(self):
        for i in range(self.rebalance_delay + 1): 
            self.x_buffer.append(self.lockin.x)
            self.y_buffer.append(self.lockin.y)
            sleep(0.5)

    #log.info('SQUID_Balancer has been called')

    def lockin_single_measure(self, lockin_amp:SR830):
        # need to add outupt conversion implementation
        x, y = lockin_amp.xy
        return x, y
    
    def lockin_buffer_measure(self, lockin_amp:SR830):

        return

    def startup(self):
        log.info('Starting dummy squid balancer')        
        self.lockin = SR830(f'GPIB{self.sr830addr}::INSTR')
        self.afg = AFG3152C(f'GPIB{self.afgaddr}::INSTR')

        self.x_buffer = FixedSizeBuffer(size=self.buffer_size)
        self.y_buffer = FixedSizeBuffer(size=self.buffer_size)

        log.info('filling the buffer')
        
        for i in range(self.x_buffer.size):
            x, y = self.lockin_single_measure(self.lockin)
            self.x_buffer.append(x)
            self.y_buffer.append(y)

        log.info('buffer filled')
        log.info('X buffer mean: {:.3g}, median: {:.3g}, std_dev: {:.3g}'.format(**self.x_buffer.compute_statistics()))
        log.info('Y buffer:')
        log.info('Y buffer mean: {:.3g}, median: {:.3g}, std_dev: {:.3g}'.format(**self.y_buffer.compute_statistics()))

        self.tek_nx = 0
        self.tek_ny = 0
        self.amp_increment = self.increment_voltage
        #self.initial_amp = self.initial_tek_amp #7.1E-3
        self.initial_amp = self.afg.ch1.amp_vrms
        self.initial_phase = self.afg.ch1.phase_rad#.value * pi/180

        #self.phasor = complex(self.tek_nx*self.amp_increment, self.tek_ny*self.amp_increment)
        self.x_comp = self.initial_amp*cos(self.initial_phase) + self.tek_nx * self.amp_increment
        self.y_comp = self.initial_amp*sin(self.initial_phase) + self.tek_ny * self.amp_increment
        self.phasor = complex(self.x_comp, self.y_comp)

        self.amp = abs(self.phasor)
        self.phase = cmath.phase(self.phasor)

        self.afg.ch1.amp_vrms = self.amp
        self.afg.ch1.phase_rad = self.phase
        log.info('The initial amp has been set to %e' %self.afg.ch1.amp_vrms)
        log.info('The initial phase has been set to %e' %self.afg.ch1.phase_deg)

        self.t0 = time.time()
    
    def data_record(self):
        x, y = self.lockin_single_measure(self.lockin)
        x_cancel = self.amp * cos(self.phase)
        y_cancel = self.amp * sin(self.phase) 
        chi_real = x + x_cancel * self.cancel_gain * self.teksr_gain
        chi_imag = y + y_cancel * self.cancel_gain * self.teksr_gain

        # self.chi_re = (lockin.x + self.cancel_gain*gain*self.amp*cos(self.phase))*self.chi_norm_constant
        # self.chi_im = (lockin.y + self.cancel_gain*gain*self.amp*sin(self.phase))*self.chi_norm_constant        
        self.data = {
            'UTC': time.time() +2082844800, #the extra constant is to match Laview's timestamp start of 01/01/1904,
            'timestamp': time.time()-self.t0,
            'X': x,
            'Y': y,
            'lockin_drive': self.lockin.sine_voltage,
            'tek_nx' : self.tek_nx,
            'tek_ny': self.tek_ny,
            'X_cancel': x_cancel,
            'Y_cancel': y_cancel,#self.phase*180/pi,
            'chi_real': chi_real,
            'chi_imag': chi_imag 
        }
        self.emit('results', self.data)        

    def adjust_cancellation(self):
        x_cancel = self.initial_amp * cos(self.initial_phase) + self.tek_nx * self.amp_increment
        y_cancel = self.initial_amp * sin(self.initial_phase) + self.tek_ny * self.amp_increment
        phasor = complex(x_cancel, y_cancel)
        amp = abs(phasor)
        phase = cmath.phase(self.phasor)
        
        self.afg.ch1.amp_vrms = amp
        self.afg.ch1.phase_rad = phase
    # def cancellation_adjust(self):

    def execute(self):
        lockin = self.lockin
        afg = self.afg
        
        while True:
            self.x_buffer.append(lockin.x)
            self.y_buffer.append(lockin.y)
            
            x_tape = np.array(self.x_buffer)
            y_tape = np.array(self.y_buffer)
            range_bound_multiplier = 1.2
            x_outof_range = np.all(np.abs(x_tape) > range_bound_multiplier * self.squid_range)
            y_outof_range = np.all(np.abs(y_tape) > range_bound_multiplier * self.squid_range)

            if x_outof_range: #x_buffer used for x balancing
                log.info('Lockin X is out of squid_range')                
                if np.median(x_tape) > 0:
                    log.info('Lockin X channel is large and positive')
                    while lockin.x > -0.8 * self.squid_range:
                        # self.x_buffer.append(lockin.x)
                        # self.y_buffer.append(lockin.y)

                        self.tek_nx = self.tek_nx + 1
                        self.adjust_cancellation()
                        self.data_record()

                        log.debug('RAISING X cancellation')
                        sleep(self.increment_sleep)
                        sleep(self.delay)
                          
                elif np.median(x_tape) < 0:
                    log.info('Lockin X channel is large and negative')
                    while lockin.x < 0.8*self.squid_range:
                        # self.x_buffer.append(lockin.x)
                        # self.y_buffer.append(lockin.y)

                        self.tek_nx = self.tek_nx-1
                        self.adjust_cancellation()
                        self.data_record()

                        log.debug('LOWERING X cancellation')
                        sleep(self.increment_sleep)
                        sleep(self.delay)   
            
                self.x_buffer.append(self.data['X'])
                self.y_buffer.append(self.data['Y'])

            if y_outof_range:
                log.info('Lockin Y channel is out of squid_range')
                #  log.info(self.y_buffer[-1*self.rebalance_delay:])
                if np.median(y_tape) > 0:
                    log.info('Lockin Y channel is large and positive')
                    while lockin.y > -0.8*self.squid_range:
                        self.x_buffer.append(lockin.x)
                        self.y_buffer.append(lockin.y)

                        self.tek_ny = self.tek_ny + 1
                        self.adjust_cancellation()
                        self.data_record()

                        log.debug('RAISING Y cancellation ')
                        sleep(self.increment_sleep)
                        sleep(self.delay)
                          
                elif np.median(y_tape) < 0:
                    log.info('Lockin Y channel is large and negative')
                    while lockin.y < 0.8*self.squid_range:
                        self.x_buffer.append(lockin.x)
                        self.y_buffer.append(lockin.y)

                        self.tek_ny = self.tek_ny-1
                        self.adjust_cancellation()
                        self.data_record()

                        log.debug('LOWERING Y cancellation')
                        sleep(self.increment_sleep)
                        sleep(self.delay)

            else:
                self.data_record()
                self.x_buffer.append(self.data['X'])
                self.y_buffer.append(self.data['Y'])  
                log.debug('Cancellation unchanged')
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
        self.setWindowTitle('SQUID dummy balancer')
        self.directory = r'D:/Data/SQUID-Dummy/'  
    
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