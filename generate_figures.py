"""
centralized figure generation file. Also acts as code documentation.
"""

# user modules
import sim_V2b as sv2

# system modules
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib



#### set up LaTeX for plots.
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'
#matplotlib.rcParams['font.size'] = 12

preamble = (r'\usepackage{amsmath}'
            r'\usepackage{siunitx}'
            r'\usepackage{bm}'
            r'\newcommand{\ve}{\varepsilon}')

matplotlib.rcParams['text.latex.preamble'] = preamble
fontsize = 12

def test():
    fig,axs = plt.subplots()
    return fig


def generate_figure(function, args, filenames, dpi=200):

    fig = function(*args)

    if type(filenames) is list:
        for name in filenames:
            fig.savefig(name,dpi=dpi)
    else:
        fig.savefig(filenames,dpi=dpi)


def main():

    # create figs directory if it doesn't exist
    if not(os.path.isdir('figs')):
        os.mkdir('figs')
    
    figures = [
        (test,[],['figs/test.pdf']),
        
    ]
    
    for fig in figures:
        generate_figure(*fig)


if __name__ == "__main__":
    main()
