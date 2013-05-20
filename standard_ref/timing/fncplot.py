# https://github.com/DropD/ethpy/blob/master/setup.py

import matplotlib.pyplot as plt

class fncplot(object):
    '''
(static class) plotting utilities for generating FNC style plots.
- xlabel: alias to pyplot.xlabel
- ylabel: places the y label on top of the y axis
in horizontal rotation using pyplot.text
- title: same as pyplot.title but shiftet up due to custom
ylabel position
'''
    @staticmethod
    def xlabel(*args, **kwargs):
        if 'offset' in kwargs:
            x += offset[0]
            y += offset[1]
        plt.xlabel(*args, **kwargs)

    @staticmethod
    def ylabel(*args, **kwargs):
        subplot = plt.gcf().get_axes()
        kwargs['transform'] = subplot[0].transAxes
        x = 0
        y = 1.04
        if 'offset' in kwargs:
            x += offset[0]
            y += offset[1]
        plt.text(x, y, *args, **kwargs)

    @staticmethod
    def title(*args, **kwargs):
        kwargs['y'] = 1.06
        plt.title(*args, **kwargs)
