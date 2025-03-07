"""
centralized figure generation file. Also acts as code documentation.
"""

import numpy as np
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


def generate_figure(function, args, filenames, dpi=200):
    # workaround for python bug where forked processes use the same random 
    # filename.
    #tempfile._name_sequence = None;

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
        (forcing_fn,[],['figs/f_forcing.pdf']),
        
    ]
    
    for fig in figures:
        generate_figure(*fig)


if __name__ == "__main__":
    main()
