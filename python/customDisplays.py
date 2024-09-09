# Some auxiliary Python functions

import sympy
sympy.init_printing()

# jupyter notebook display library
from IPython.display import display as IPyDisplay
from IPython.display import display, Latex

PyPrint = print


def doNothing(*ka, **kwa): pass


## Decomment the following lines if you want to disable printing:
## useful if the rendering is very slow.
# print = doNothing
# display = doNothing

## To restore the original printing functions, do this:
# print = PyPrint
# display = IPyDisplay

# Custom printing functions

def latexColor(s, color):
    return "{\\color{%s}{" % color + s + '}}'


def vPhantom(obj):
    objt = type(obj)
    if objt in (list, tuple, set, dict):
        res = ''
        objlist = list(obj.keys()) if objt == dict else tuple(obj)
        for o in objlist:
            res += " " + vPhantom(o)
            if objt == dict:
                res += " " + vPhantom(obj[o])
    else:
        res = "\\vphantom{%s}" % sympy.latex(obj)
    return res


def toCustomLatex(
        obj,
        recursive=False,
        delimiters=True,
        separator=',',
        spacer='\\quad',
        color='black'
):
    objt = type(obj)

    Color = lambda x: latexColor(x, color)
    text2latex = lambda x: Color("\\text{%s}" % str(x))

    if type(obj) == str:
        return text2latex(obj)

    objlist = list(obj.keys()) if objt == dict else tuple(obj)
    latex_obj = ''

    toRecLatex = lambda x: toCustomLatex(
        x,
        recursive=recursive,
        delimiters=True,
        color=color
    )

    for ii in range(len(objlist)):

        if type(objlist[ii]) in (set, dict, list, tuple) and recursive:
            toLatex = toRecLatex
        elif type(objlist[ii]) == str:
            toLatex = text2latex
        else:
            toLatex = lambda x: Color(sympy.latex(x))

        try:
            o = toLatex(objlist[ii])
        except:
            o = text2latex(o)

        if objt == dict:
            try:
                o += ": " + toLatex(obj[objlist[ii]])
            except:
                o += ": " + text2latex(obj[objlist[ii]])

        sep = (ii + 1 < len(objlist)) * (Color(separator) + spacer + '\n')
        latex_obj += o + sep

    if delimiters:
        ldelims, rdelims = {
            set: "\{ \}",
            dict: "\{ \}",
            list: "[ ]",
            tuple: "( )"
        }[objt].split(" ")

        ldelim = \
            "\\left " + ldelims + vPhantom(obj) + "\\right .\n"
        rdelim = \
            "\\left . " + vPhantom(obj) + "\\right " + rdelims + "\n"

        latex_obj = Color(ldelim) + latex_obj + Color(rdelim)

    return latex_obj


def customDisplay(
        obj,
        recursive=False,
        color='black',
        textcolor='brown',
        **args
):
    objt = type(obj)

    Color = lambda x: latexColor(x, color)

    if objt in (set, dict, list, tuple):
        try:
            latex_obj = toCustomLatex(
                obj,
                recursive=recursive,
                color=color,
                **args
            )
        except:
            latex_obj = toCustomLatex(str(obj))

        latex_obj = "\\begin{align*}\n\\begin{autobreak}\n" \
                    + latex_obj \
                    + "\\end{autobreak}\n\\end{align*}\n" \

        display(Latex(latex_obj))
        return

    if objt == str:
        display(Latex("$%s$" % latexColor("\\text{%s}" % obj, textcolor)))
        return

    try:
        obj_latex = sympy.latex(obj)
        display(Latex("$$%s$$" % Color(obj_latex)))
    except:
        try:
            customDisplay(str(obj))
        except:
            display(obj)


def blockDisplay(*args, **kwargs):
    for a in args:
        customDisplay(a, **kwargs)


# Precision used for the definition of numeric values in Sympy
# equations.
precision = 50


def NVal(val):
    """
        Shortcut to create Sympy values with arbitrary precision.
    """
    return sympy.N(val, precision)


def chop(obj, low_precision=2):
    """
        Truncate the very long floating point values of the Sympy
        object obj. This function is meant to be used only for
        displaying, to make Sympy objects readable.
    """
    if isinstance(obj, dict):
        return {kk: chop(obj[kk]) for kk in obj.keys()}
    for typ in [list, tuple, set]:
        if isinstance(obj, typ):
            return typ(chop(o) for o in obj)
    if type(obj) in (sympy.Float, float):
        exponent = -sympy.floor(sympy.log(sympy.Abs(obj)) \
                                / sympy.log(10.)) + low_precision
        return sympy.N(sympy.Integer(obj * 10 ** exponent) \
                       / 10 ** exponent, precision)
    sfloat = sympy.Wild(
        'sfloat',
        properties=[lambda x: isinstance(x, sympy.Float)]
    )
    if str(type(obj))[8:13] == 'sympy':
        return sympy.N(obj.replace(
            sfloat,
            lambda sfloat: chop(sfloat))
        )
    return obj
