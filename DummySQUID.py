import logging
import numpy as np
from utils import FixedSizeBuffer

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
from pymeasure.display.widgets import PlotWidget, LogWidget #, PlotWidget_datetimeaxis
from pymeasure.experiment import Procedure, Results, unique_filename
from pymeasure.experiment import IntegerParameter, FloatParameter, Parameter

from PyQt5 import QtWidgets
from pyqtgraph import DateAxisItem

import pyvisa

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

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
    increment_sleep = FloatParameter('Increment sleep', units = 's', default = 0.1)
    # drive_V = FloatParameter('drive voltage', units = 'V', default = 3.25)
    cancel_gain = FloatParameter('cancellation gain', default = 1)
    squid_range = FloatParameter('SQUID output/Lockin input limit',units= 'V', default = 10E-3) #the max range on the self.lockin
    increment_voltage = FloatParameter('Increment voltage', units = 'V', default = 1E-3)
    buffer_size = IntegerParameter('Re-balance buffer count', default = 30)
    


    sr830addr = Parameter('SR830 address', default='2::9')
    afgaddr = Parameter('AFG address', default='2::11')

    teksr_gain = FloatParameter('Tek/lock-in gain', default = 1)

    inputs_displays = [ 
                    'delay',
                    'increment_sleep',
                    # 'drive_V',
                    'cancel_gain',
                    'increment_voltage',
                    'squid_range',
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
        'signal_x', 'signal_y'
        ] 
                
    # def find_gain(self):

    #     gain = 1  #just some code to not cause an error 
        
    #     # zeros the input on the lock-in x channel with the tektronix and tells you what the gain is
    #     # gain = self.lockin.x/tek_amp
    #     self.amp = self.initial_amp
    #     self.afg.ch1.amp_vrms = self.amp
    #     self.afg.ch1.disable()  #set the amplitude to the lowest setting and disable to counterbalance
    #     signal = self.lockin.x  
        
    #     self.afg.ch1.enable()
        
    #     while abs(self.lockin.x) > 0.01*self.squid_range:
    #         tek_nx = tek_nx + 1
    #         phasor = complex(tek_nx * self.increment_voltage, tek_ny * self.increment_voltage)
    #         amp = self.initial_amp + abs(phasor)
    #         phase = self.initial_phase + cmath.phase(phasor)
            
    #         self.afg.ch1.amp_vrms = self.amp; #print(afg.ch1.amp_vrms)
    #         self.afg.ch1.phase_rad = self.phase
    
    def SR830_single_measure(self, lockin_amp:SR830):
        # need to add outupt conversion implementation
        x, y = lockin_amp.xy
        return x, y
    
    def lockin_cancelbuffer_measure(self, lockin_amp:SR830):

        return

    def tek_signal_query(self, tek_afg:AFG3152C, phase_units='rad'):
        # need to include assertion error for phase_units
        state = True
        while state:
            try: 
                amp = tek_afg.ch1.amp_vrms
                if phase_units=='rad':
                    phase = tek_afg.ch1.phase_rad
                else:
                    phase = tek_afg.ch1.phase_deg
                state = False
            except Exception as e:
                log.warning(e)
                time.sleep(0.1)
        return amp, phase

    def data_record(self):
        x, y = self.SR830_single_measure(self.lockin)
        # amp = self.afg.ch1.amp_vrms
        # phase = self.afg.ch1.phase_rad
        amp, phase = self.tek_signal_query(self.afg)
        x_cancel = amp * cos(phase)
        y_cancel = amp * sin(phase) 
        signal_x = x + x_cancel * self.cancel_gain * self.teksr_gain
        signal_y = y + y_cancel * self.cancel_gain * self.teksr_gain

       
        self.data = {
            'UTC': time.time() + 2082844800, #the extra constant is to match Laview's timestamp start of 01/01/1904,
            'timestamp': time.time()-self.t0,
            'X': x,
            'Y': y,
            'lockin_drive': self.lockin.sine_voltage,
            'tek_nx' : self.tek_nx,
            'tek_ny': self.tek_ny,
            'X_cancel': x_cancel,
            'Y_cancel': y_cancel,#self.phase*180/pi,
            'signal_x': signal_x,
            'signal_y': signal_y 
        }
        self.emit('results', self.data)        

    def adjust_cancellation(self):
        x_cancel = self.initial_amp * cos(self.initial_phase) + self.tek_nx * self.increment_voltage
        y_cancel = self.initial_amp * sin(self.initial_phase) + self.tek_ny * self.increment_voltage
        phasor = complex(x_cancel, y_cancel)
        amp = abs(phasor)
        phase = cmath.phase(phasor)
        self.afg.ch1.amp_vrms = amp
        self.afg.ch1.phase_rad = phase

    def log_cancellation(self):
        Xlog, Ylog = self.data['X_cancel'], self.data['Y_cancel']
        Rlog = (Xlog**2 + Ylog**2)**0.5
        philog = np.arctan(Ylog/Xlog) * 180 / pi

        # `log` is a GLOBAL variable
        log.info('Cancellation (R, phi) = ({:.3e}V, {:.3g}deg)'.format(Rlog, philog))
        log.info('Cancellation (X, Y) = ({:.3e}V, {:.3e}V)'.format(Xlog, Ylog))        
        return        

    def startup(self):
        log.info('Starting dummy squid balancer')        
        self.lockin = SR830(f'GPIB{self.sr830addr}::INSTR')
        self.afg = AFG3152C(f'GPIB{self.afgaddr}::INSTR')
        log.info(self.lockin.id); log.info(self.afg.id)

        self.x_cancelbuffer = FixedSizeBuffer(size=self.buffer_size)
        self.y_cancelbuffer = FixedSizeBuffer(size=self.buffer_size)

        log.info('filling the buffer')
        
        for _ in range(self.buffer_size):
            x, y = self.SR830_single_measure(self.lockin)
            # log.info(x); log.info(y)
            self.x_cancelbuffer.append(x)
            self.y_cancelbuffer.append(y)

        log.info('buffer filled')

        # self.initial_amp = self.afg.ch1.amp_vrms
        # self.initial_phase = self.afg.ch1.phase_rad
        self.initial_amp, self.initial_phase = self.tek_signal_query(self.afg, 'deg')
        xcancel_initial = self.initial_amp * cos(self.initial_phase)
        ycancel_initial = self.initial_amp * sin(self.initial_phase)
        self.tek_nx = 0
        self.tek_ny = 0
        log.info('Initial (R, phi): ({:}V, {:.3g}deg)'.format(self.initial_amp, self.initial_phase))
        log.info('Initial (X, Y): ({:.3g}V, {:.3g}V)'.format(xcancel_initial, ycancel_initial))
        self.adjust_cancellation()

        self.t0 = time.time()

    def execute(self):
                
        while True:
            self.x_cancelbuffer.append(self.lockin.x)
            self.y_cancelbuffer.append(self.lockin.y)
            
            x_tape = self.x_cancelbuffer.get_buffer()
            y_tape = self.y_cancelbuffer.get_buffer()
            range_bound_multiplier = 1.2
            x_outof_range = np.all(np.abs(x_tape) > range_bound_multiplier * self.squid_range)
            y_outof_range = np.all(np.abs(y_tape) > range_bound_multiplier * self.squid_range)

            if x_outof_range: #x_cancelbuffer used for x balancing
                log.info('Lockin X is out of squid_range')                
                if np.median(x_tape) > 0:
                    log.info('Lockin X channel is large and positive')
                    while self.lockin.x > -0.8 * self.squid_range:
                        self.tek_nx += 1

                        try:
                            self.adjust_cancellation()
                        except Exception as e:
                            log.warning(f"Exception occurred during adjustment: {e}. Setting amplitude to floor limit.")
                            self.afg.ch1.amp_vrms = 7e-3  # Set to floor limit
                            continue  # Continue the loop
                        
                        self.data_record()

                        log.debug('RAISING X cancellation')
                        sleep(self.increment_sleep)
                        sleep(self.delay)
                        if self.should_stop():
                            log.warning('Caught the stop flag in the procedure')
                            break

                    self.log_cancellation()
                    # log.info('Cancellation (R, phi) = ({:.3e}V, {:.3g}deg)'.format(*self.tek_signal_query(self.afg, 'deg')))
                    # log.info('Cancellation (X, Y) = ({:.3e}V, {:.3e}V)'.format(self.data['X_cancel'], self.data['Y_cancel']))

                elif np.median(x_tape) < 0:
                    log.info('Lockin X channel is large and negative')
                    while self.lockin.x < 0.8*self.squid_range:

                        self.tek_nx -= 1
                        try:
                            self.adjust_cancellation()
                        except Exception as e:
                            log.warning(f"Exception occurred during adjustment: {e}. Setting amplitude to floor limit.")
                            self.afg.ch1.amp_vrms = 7e-3  # Set to floor limit
                            continue  # Continue the loop                  
                        self.data_record()

                        log.debug('LOWERING X cancellation')
                        sleep(self.increment_sleep)
                        sleep(self.delay)   
                        if self.should_stop():
                            log.warning('Caught the stop flag in the procedure')
                            break

                    self.log_cancellation()
                    # log.info('Cancellation (R, phi) = ({:.3e}V, {:.3g}deg)'.format(*self.tek_signal_query(self.afg, 'deg')))
                    # log.info('Cancellation (X, Y) = ({:.3e}V, {:.3e}V)'.format(self.data['X_cancel'], self.data['Y_cancel']))

                self.x_cancelbuffer.append(self.data['X'])
                self.y_cancelbuffer.append(self.data['Y'])

            if y_outof_range:
                log.info('Lockin Y channel is out of squid_range')
                if np.median(y_tape) > 0:
                    log.info('Lockin Y channel is large and positive')
                    while self.lockin.y > -0.8*self.squid_range:
                        self.x_cancelbuffer.append(self.lockin.x)
                        self.y_cancelbuffer.append(self.lockin.y)

                        self.tek_ny += 1
                        try:
                            self.adjust_cancellation()
                        except Exception as e:
                            log.warning(f"Exception occurred during adjustment: {e}. Setting amplitude to floor limit.")
                            self.afg.ch1.amp_vrms = 7e-3  # Set to floor limit
                            continue  # Continue the loop               
                        self.data_record()

                        log.debug('RAISING Y cancellation ')
                        sleep(self.increment_sleep)
                        sleep(self.delay)
                        if self.should_stop():
                            log.warning('Caught the stop flag in the procedure')
                            break

                    self.log_cancellation()
                    # log.info('Cancellation (R, phi) = ({:.3e}V, {:.3g}deg)'.format(*self.tek_signal_query(self.afg, 'deg')))
                    # log.info('Cancellation (X, Y) = ({:.3e}V, {:.3e}V)'.format(self.data['X_cancel'], self.data['Y_cancel']))

                elif np.median(y_tape) < 0:
                    log.info('Lockin Y channel is large and negative')
                    while self.lockin.y < 0.8*self.squid_range:
                        self.x_cancelbuffer.append(self.lockin.x)
                        self.y_cancelbuffer.append(self.lockin.y)

                        self.tek_ny -= 1
                        try:
                            self.adjust_cancellation()
                        except Exception as e:
                            log.warning(f"Exception occurred during adjustment: {e}. Setting amplitude to floor limit.")
                            self.afg.ch1.amp_vrms = 7e-3  # Set to floor limit
                            continue  # Continue the loop                        self.data_record()

                        log.debug('LOWERING Y cancellation')
                        sleep(self.increment_sleep)
                        sleep(self.delay)
                        if self.should_stop():
                            log.warning('Caught the stop flag in the procedure')
                            break
                    
                    self.log_cancellation()
                    # log.info('Cancellation (R, phi) = ({:.3e}V, {:.3g}deg)'.format(*self.tek_signal_query(self.afg, 'deg')))
                    # log.info('Cancellation (X, Y) = ({:.3e}V, {:.3e}V)'.format(self.data['X_cancel'], self.data['Y_cancel']))

            else:
                self.data_record()
                self.x_cancelbuffer.append(self.data['X'])
                self.y_cancelbuffer.append(self.data['Y'])  
                # log.debug('Cancellation unchanged')
                sleep(self.delay)
            
            if self.should_stop():
                log.warning('Caught the stop flag in the procedure')
                break
 
class lockinOutput(ManagedWindow): 
    def __init__(self):
        # ChiPlot = PlotWidget( #PlotWidget_datetimeaxis(
        #     name = 'Chi', 
        #     columns = SQUID_Balancer_v2.DATA_COLUMNS,
        #     x_axis = 'UTC',
        #     y_axis = 'chi\'')

        #ChiPlot.plot.setAxisItems({'bottom': DateAxisItem(utcOffset= 2082844800)})

        # TempPlot = PlotWidget_datetimeaxis(
        #     name = 'Temperature data', 
        #     columns = SQUID_Balancer_v2.DATA_COLUMNS, 
        #     x_axis = 'UTC', 
        #     y_axis = 'Temperature')

        super().__init__(
            procedure_class=SQUID_Balancer_v2,
            inputs=  SQUID_Balancer_v2.inputs_displays,
            displays = SQUID_Balancer_v2.inputs_displays,
            x_axis = 'UTC',
            y_axis = 'X',
            # directory_input= True,
            # widget_list = (
            #     # TempPlot, 
            #     ChiPlot,
            # )
        )
        self.setWindowTitle('SQUID dummy balancer')
        self.directory = r'D:/Data/SQUID-Dummy/'  
    
    def queue(self):
        directory = self.directory
        filename = unique_filename(directory, prefix = 'Dummy')
        procedure = self.make_procedure()
        results = Results(procedure, filename)
        experiment = self.new_experiment(results)
        self.manager.queue(experiment)
    


if __name__ == '__main__':
    
    rm = pyvisa.ResourceManager()
    gpib_list = rm.list_resources()
    print(pymeasure.__version__)

    # app = QtGui.QApplication(sys.argv)
    app = QtWidgets.QApplication([])
    window = lockinOutput()
    window.show()
    sys.exit(app.exec_())