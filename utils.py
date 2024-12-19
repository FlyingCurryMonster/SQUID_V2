import numpy as np
from collections import deque

class FixedSizeBuffer:
    def __init__(self, size=30):
        self.buffer = deque(maxlen=size)

    def append(self, element):
        self.buffer.append(element)

    def compute_statistics(self):
        data = np.array(self.buffer)
        return {
            "mean": np.mean(data),
            "median": np.median(data),
            "std_dev": np.std(data),
        }


