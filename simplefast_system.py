from compartment_manager import SimpleFastCompartmentManager
import numpy as np

time = 400

if __name__ == "__main__":
    sfcm = SimpleFastCompartmentManager((19, 1, 0, 6), 18, 1)
    sfcm.simComparts(time)


