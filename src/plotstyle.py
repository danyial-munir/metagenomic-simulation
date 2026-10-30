"""One figure style for the whole project; Okabe-Ito colour-blind-safe palette."""
import matplotlib as mpl

PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#000000"]


def set_style():
    mpl.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.prop_cycle": mpl.cycler(color=PALETTE),
        "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    })
