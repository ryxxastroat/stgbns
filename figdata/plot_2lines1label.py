import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple

class HandlerTupleVertical(HandlerTuple):
    def __init__(self, **kwargs):
        HandlerTuple.__init__(self, **kwargs)

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        # How many lines are there.
        numlines = len(orig_handle)
        handler_map = legend.get_legend_handler_map()

        # divide the vertical space where the lines will go
        # into equal parts based on the number of lines
        height_y = (height / numlines)

        leglines = []
        for i, handle in enumerate(orig_handle):
            handler = legend.get_legend_handler(handler_map, handle)

            legline = handler.create_artists(legend, handle,
                                             xdescent,
                                             (2*i + 1)*height_y,
                                             width,
                                             2*height,
                                             fontsize, trans)
            leglines.extend(legline)

        return leglines
        
        
        
fig, ax1 = plt.subplots(1, 1)

# two legend keys for a single entry
p2, = ax1.plot([3, 4], [2, 3], '-')
p1, = ax1.plot([1, 2], [5, 6], '--')
# `plot` returns a list, but we want the handle - thus the comma on the left

# Assign two of the handles to the same legend entry by putting them in a tuple
# and using a generic handler map (which would be used for any additional
# tuples of handles like (p1, p3)).
l = ax1.legend([(p1, p2)], ['some label'],
               handler_map={tuple: HandlerTupleVertical()})
               
               
plt.show()               
