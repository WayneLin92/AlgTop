"""This modules provide some functions to be used in Jupyter notebook."""
import IPython.display as display
from inspect import signature


def help_html(fn):
    result = f"{fn.__name__}{str(signature(fn))}<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;{fn.__doc__}"
    return display.HTML(result)
