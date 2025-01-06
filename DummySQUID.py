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
    
    delay = FloatParameter('delay', units='s', default=1)
    increment_sleep = FloatParameter('Increment sleep', units = 's', default=0.1)
    
    buffer_size = IntegerParameter('Re-balance buffer count', default=30)
    squid_range = FloatParameter('SQUID output/Lockin input limit',units='V', default=10E-3) #the max range on the self.lockin
    increment_voltage = FloatParameter('Increment voltage', units='V', default=1E-3)
    
    cancel_gain = FloatParameter('cancellation gain', default=1)
    teksr_gain = FloatParameter('Tek/lock-in gain', default=2)

    sr830addr = Parameter('SR830 address', default='2::9')
    afgaddr = Parameter('AFG address', default='2::11')

    inputs_displays = [ 
        'delay',
        'increment_sleep',
        'buffer_size',
        'squid_range',
        'increment_voltage',
        'cancel_gain',
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

    def SR830_single_measure(self, lockin_amp:SR830):
        # need to add outupt conversion implementation
        x, y = lockin_amp.xy
        return x, y

    def tek_signal_query(self, tek_afg:AFG3152C, phase_units='rad'):
        # need to include assertion error for phase_units
        assert phase_units in ['rad', 'deg']

        state = True
        while state:
            try: 
                amp = tek_afg.ch1.amp_vrms
                if phase_units=='rad':
                    phase = tek_afg.ch1.phase_rad
                else:
                    phase = tek_afg.ch1.phase_deg
                
                assert isinstance(amp, float), f"Amplitude is not a float: {type(amp)}"
                assert isinstance(phase, float), f"Phase is not a float: {type(phase)}"

                state = False

            except AssertionError as ae:
                log.warning(f"Assertion failed: {ae}. Retrying...")
                time.sleep(0.1)  
              
            except Exception as e:
                log.warning(e)
                time.sleep(0.1)
        return amp, phase

    def data_record(self):
        x, y = self.SR830_single_measure(self.lockin)

        # not sure if I should update the buffer after a measurement
        # self.x_cancelbuffer.append(x)
        # self.y_cancelbuffer.append(y)
        
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

    def cancellation_incrementer(self, lockin_amp: SR830, component, direction):
        """
        Incrementally adjusts cancellation for a specific component and direction.

        Args:
            lockin_amp (SR830): The lock-in amplifier object.
            component (int): The component to adjust (0 for X, 1 for Y).
            direction (int): The direction of adjustment (1 for positive, -1 for negative).
        """
        try:
            # Validate inputs
            assert component in [0, 1], f"Invalid component: {component}. Must be 0 (X) or 1 (Y)."
            assert direction in [1, -1], f"Invalid direction: {direction}. Must be 1 or -1."
        except AssertionError as ae:
            log.error(f"Input validation failed: {ae}")
            return  # Abort if inputs are invalid
        
        action = "RAISED" if direction == 1 else "LOWERED"
        component_name = "X" if component == 0 else "Y"   
        direction_name = "INCREASING" if direction == 1 else "DECREASING"

        log.info(f"Adjusting cancellation for component {component_name} in {direction_name} direction.")
        # tape_buffer = self.x_cancelbuffer if component == 0 else self.y_cancelbuffer
        squid_range_limit = -0.8 * self.squid_range if direction == 1 else 0.8 * self.squid_range

        while (lockin_amp.xy[component] > squid_range_limit) if direction == 1 else (lockin_amp.xy[component] < squid_range_limit):
            if component == 0:  # X component
                self.tek_nx += direction
            else:  # Y component
                self.tek_ny += direction
            
            try:
                self.adjust_cancellation()
            except Exception as e:
                log.warning(f"Error during cancellation adjustment: {e}. Retrying...")
                continue

            sleep(self.increment_sleep)
            self.data_record()
            sleep(self.delay)

            if self.should_stop():
                log.warning('Caught the stop flag in the procedure')
                break

        log.info(f"{action} {component_name} cancellation")
        
        Xlog, Ylog = self.data['X_cancel'], self.data['Y_cancel']
        Rlog = (Xlog**2 + Ylog**2)**0.5
        philog = np.arctan(Ylog/Xlog) * 180 / pi

        # `log` is a GLOBAL variable
        log.info('Cancellation (R, phi) = ({:.3e}V, {:.3g}deg)'.format(Rlog, philog))
        log.info('Cancellation (X, Y) = ({:.3e}V, {:.3e}V)'.format(Xlog, Ylog))      

    def adjust_cancellation(self):
        x_cancel = self.initial_amp * cos(self.initial_phase) + self.tek_nx * self.increment_voltage
        y_cancel = self.initial_amp * sin(self.initial_phase) + self.tek_ny * self.increment_voltage
        phasor = complex(x_cancel, y_cancel)
        amp = abs(phasor)
        phase = cmath.phase(phasor)
        
        state = True
        while state:
            try: 
                self.afg.ch1.amp_vrms = amp
                self.afg.ch1.phase_rad = phase
                state = False
            except Exception as e:
                log.warning(f"Exception occurred during adjustment: {e}. Setting amplitude to floor limit.")
                self.afg.ch1.amp_vrms = 7e-3
                continue   

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
        sleep(self.delay * 3)
        while True:
            self.x_cancelbuffer.append(self.lockin.x)
            self.y_cancelbuffer.append(self.lockin.y)

            x_tape = self.x_cancelbuffer.get_buffer()
            y_tape = self.y_cancelbuffer.get_buffer()
            range_bound_multiplier = 1.2
            x_outof_range = np.all(np.abs(x_tape) > range_bound_multiplier * self.squid_range)
            y_outof_range = np.all(np.abs(y_tape) > range_bound_multiplier * self.squid_range)

            if x_outof_range:
                log.info('Lockin X is out of squid_range')
                if np.median(x_tape) > 0:
                    self.cancellation_incrementer(self.lockin, 0, 1)  # Increase X cancellation
                elif np.median(x_tape) < 0:
                    self.cancellation_incrementer(self.lockin, 0, -1)  # Decrease X cancellation

            if y_outof_range:
                log.info('Lockin Y is out of squid_range')
                if np.median(y_tape) > 0:
                    self.cancellation_incrementer(self.lockin, 1, 1)  # Increase Y cancellation
                elif np.median(y_tape) < 0:
                    self.cancellation_incrementer(self.lockin, 1, -1)  # Decrease Y cancellation

            else:
                self.data_record()
                self.x_cancelbuffer.append(self.data['X'])
                self.y_cancelbuffer.append(self.data['Y'])
                sleep(self.delay)

            if self.should_stop():
                log.warning('Caught the stop flag in the procedure')
                break

    # def execute(self):
    #     sleep(self.delay * 3)
    #     while True:
    #         self.x_cancelbuffer.append(self.lockin.x)
    #         self.y_cancelbuffer.append(self.lockin.y)
            
    #         x_tape = self.x_cancelbuffer.get_buffer()
    #         y_tape = self.y_cancelbuffer.get_buffer()
    #         range_bound_multiplier = 1.2
    #         x_outof_range = np.all(np.abs(x_tape) > range_bound_multiplier * self.squid_range)
    #         y_outof_range = np.all(np.abs(y_tape) > range_bound_multiplier * self.squid_range)

    #         if x_outof_range: #x_cancelbuffer used for x balancing
    #             log.info('Lockin X is out of squid_range')                
    #             if np.median(x_tape) > 0:
    #                 log.info('Lockin X channel is large and positive')
    #                 while self.lockin.x > -0.8 * self.squid_range:
    #                     self.tek_nx += 1

    #                     self.adjust_cancellation()
    #                     sleep(self.increment_sleep)
    #                     self.data_record()
    #                     sleep(self.delay)
    #                     if self.should_stop():
    #                         log.warning('Caught the stop flag in the procedure')
    #                         break

    #                 log.info('RAISED X cancellation')
    #                 self.log_cancellation()

    #             elif np.median(x_tape) < 0:
    #                 log.info('Lockin X channel is large and negative')
    #                 while self.lockin.x < 0.8*self.squid_range:

    #                     self.tek_nx -= 1
    #                     self.adjust_cancellation()
    #                     sleep(self.increment_sleep)
    #                     self.data_record()
    #                     sleep(self.delay)

    #                     if self.should_stop():
    #                         log.warning('Caught the stop flag in the procedure')
    #                         break
                    
    #                 log.info('LOWERED X cancellation')
    #                 self.log_cancellation()

    #             self.x_cancelbuffer.append(self.data['X'])
    #             self.y_cancelbuffer.append(self.data['Y'])

    #         if y_outof_range:
    #             log.info('Lockin Y channel is out of squid_range')
    #             if np.median(y_tape) > 0:
    #                 log.info('Lockin Y channel is large and positive')
    #                 while self.lockin.y > -0.8*self.squid_range:
    #                     self.x_cancelbuffer.append(self.lockin.x)
    #                     self.y_cancelbuffer.append(self.lockin.y)

    #                     self.tek_ny += 1
    #                     self.adjust_cancellation()
    #                     sleep(self.increment_sleep)
    #                     self.data_record()
    #                     sleep(self.delay)
                        
    #                     if self.should_stop():
    #                         log.warning('Caught the stop flag in the procedure')
    #                         break
                    
    #                 log.info('RAISED Y cancellation ')
    #                 self.log_cancellation()

    #             elif np.median(y_tape) < 0:
    #                 log.info('Lockin Y channel is large and negative')
    #                 while self.lockin.y < 0.8*self.squid_range:
    #                     self.x_cancelbuffer.append(self.lockin.x)
    #                     self.y_cancelbuffer.append(self.lockin.y)

    #                     self.tek_ny -= 1
    #                     self.adjust_cancellation()
    #                     sleep(self.increment_sleep)
    #                     self.data_record()

    #                     sleep(self.delay)
    #                     if self.should_stop():
    #                         log.warning('Caught the stop flag in the procedure')
    #                         break

    #                 log.info('LOWERED Y cancellation')                    
    #                 self.log_cancellation()

    #             self.x_cancelbuffer.append(self.data['X'])
    #             self.y_cancelbuffer.append(self.data['Y'])

    #         else:
    #             self.data_record()
    #             self.x_cancelbuffer.append(self.data['X'])
    #             self.y_cancelbuffer.append(self.data['Y'])  
    #             sleep(self.delay)
            
    #         if self.should_stop():
    #             log.warning('Caught the stop flag in the procedure')
    #             break
 
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
        self.file_input.filename_fixed = False
    
    def queue(self):
        directory = self.directory
        filename = unique_filename(directory, prefix=f'{self.filename}_')
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