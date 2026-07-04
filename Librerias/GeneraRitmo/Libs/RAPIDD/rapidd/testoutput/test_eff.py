import ctypes
import numpy as np
import os

from Final.Librerias.GeneraRitmo.Libs.RAPIDD.rapidd.core import _crapidd
from Final.Librerias.GeneraRitmo.Libs.RAPIDD.rapidd.core import base_dir
from scipy.interpolate import interp1d

from Final.Librerias.GeneraRitmo.Libs.RAPIDD.rapidd.experiments import read_efficiency, efficiency_fn


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    