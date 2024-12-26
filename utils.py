import numpy as np

class FixedSizeBuffer:
    def __init__(self, size):
        self.size = size
        self.buffer = np.empty(size, dtype=float)  # Pre-allocate the buffer
        self.index = 0
        self.full = False

    def append(self, element):
        """Add a new element to the buffer."""
        self.buffer[self.index] = element
        self.index = (self.index + 1) % self.size  # Wrap around when the buffer is full
        if self.index == 0:
            self.full = True

    def get_buffer(self):
        """Return the current buffer as a NumPy array."""
        if self.full:
            return self.buffer
        return self.buffer[:self.index]  # Only return valid data if buffer isn't full




# import numpy as np
# from collections import deque

# class FixedSizeBuffer:
#     def __init__(self, size=30):
#         self.buffer = deque(maxlen=size)
#         self.size = size

#     def append(self, element):
#         self.buffer.append(element)

#     def compute_statistics(self):
#         data = np.array(self.buffer)
#         return {
#             "mean": np.mean(data),
#             "median": np.median(data),
#             "std_dev": np.std(data),
#         }

#     # def compare_to_number(self, number, comparison_target="all"):
#     #     """Compare a number to elements, mean, median, or std_dev of the buffer.
        
#     #     Args:
#     #         number (float): The number to compare against.
#     #         comparison_target (str): What to compare the number to. Options are:
#     #             - "all": Compare the number to all elements in the buffer.
#     #             - "mean": Compare the number to the mean of the buffer.
#     #             - "median": Compare the number to the median of the buffer.
#     #             - "std_dev": Compare the number to the standard deviation of the buffer.
        
#     #     Returns:
#     #         bool: True if the number satisfies the condition, False otherwise.
#     #     """
#     #     data = np.array(self.buffer)
#     #     stats = self.compute_statistics()

#     #     if comparison_target == "all":
#     #         return all(data == number)  # Adjust for your specific condition (e.g., >, <, etc.)
#     #     elif comparison_target == "mean":
#     #         return number == stats["mean"]  # Adjust for your specific condition (e.g., >, <, etc.)
#     #     elif comparison_target == "median":
#     #         return number == stats["median"]  # Adjust for your specific condition (e.g., >, <, etc.)
#     #     elif comparison_target == "std_dev":
#     #         return number == stats["std_dev"]  # Adjust for your specific condition (e.g., >, <, etc.)
#     #     else:
#     #         raise ValueError("Invalid comparison_target. Choose from 'all', 'mean', 'median', or 'std_dev'.")
