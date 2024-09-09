import sympy

##  Substitution functions

def sortSubsDict(subs_dict):
    """
        Given a dictionary of unsorted substitution dictionaries,
        return a sorted list of substitution dictionaries without
        conflicts.

        A conflict happens when one of the terms to be replaced is
        contained in another term to be replaced.
        For example, the expression 1 + G(t + 2) contains the
        expression t + 2. If both must be replaced at the same
        time, there is a conflict, since replacing t + 2 may
        make impossible to replace 1 + G(t + 2).

        In this case, the substitutions involving "bigger terms"
        must be applied first.

        The returned list of dictionaries is built to guarantee
        that there is no conflict between substitutions of the
        same dictionary, and the expressions involved in the
        substitutions of each dictionary are not contained
        in the subsequent dictionaries of the returned list.

        Accepted types: dictionary
    """
    dict_list = []
    keys = list(subs_dict.keys())
    while keys:
        dict_list += [dict()]
        for s in keys:
            if sum(len(x.find(s)) for x in keys) == 1:
                dict_list[-1].update({s: subs_dict[s]})
        for s in dict_list[-1].keys():
            keys.remove(s)
    return dict_list


def subs_aux(
        obj,
        subs_dict_or_list,
        RHSonly=False,
        simplify=True,
        doit=True,
        recursive=False,
        maxRecursionDepth=10
):
    """
        Apply the given list of substitutions subs_dict_or_list
        to the object obj.

        Accepted types for the object:
            sympy equations or expressions or anything
                having the .subs() method,
            dict, list, tuple, set

        Accepted types for the substitutions:
            list of pairs,
            dictionary of substitutions

        An object of the same type of obj is returned.
        If the object is of type dictionary, the substitutions
        are NOT applied to its keys.

        If RHSonly is True and the object is of type sympy.Eq,
        the substitutions are applied only to the right hand
        side of the equation.
        In any case, simplify() and doit() are applied separately
        to each side of an equation.

        If an equation becomes trivially True, it is removed from
        the dictionary or list or set or tuple.
    """

    if type(subs_dict_or_list) == dict:
        subs_list = [(kk, subs_dict_or_list[kk])
                     for kk in subs_dict_or_list.keys()]
    else:
        subs_list = list(subs_dict_or_list)

    subss_aux = lambda x: subs_aux(
        x,
        subs_list,
        RHSonly=RHSonly,
        simplify=simplify,
        doit=doit,
        recursive=recursive
    )

    notTrivial = \
        lambda x: not isinstance(x, sympy.boolalg.BooleanTrue)

    objt = type(obj)
    if objt == dict:
        res = dict()
        for kk in obj.keys():
            objsub = subss_aux(obj[kk])
            if notTrivial(objsub):
                res.update({kk: objsub})
    elif objt in (list, tuple, set):
        res = objt()
        for kk in obj:
            skk = subss_aux(kk)
            if notTrivial(skk):
                if objt == list:
                    res.append(skk)
                elif objt == tuple:
                    res += (skk,)
                else:
                    res.add(skk)
    else:
        if isinstance(objt, sympy.Eq):
            if RHSonly:
                lhs = obj.args[0]
            else:
                lhs = subss_aux(obj.args[0])
            rhs = subss_aux(obj.args[1])
            res = sympy.Eq(lhs, rhs)
        else:
            res = obj.subs(subs_list, simplify=False)
            if recursive:
                oldres = obj
                for ii in range(maxRecursionDepth):
                    if res != oldres:
                        oldres = res
                        res = res.subs(subs_list, simplify=False)
                    else:
                        break
            if simplify:
                res = res.simplify()
            if doit:
                res = res.doit()
    return res


def subs(
        obj,
        subs_dict,
        RHSonly=False,
        simplify=True,
        doit=True,
        recursive=False,
        maxRecursionDepth=10
):
    res = obj
    subs_dict_list = sortSubsDict(subs_dict)
    for sd in subs_dict_list:
        res = subs_aux(
            res,
            sd,
            RHSonly=RHSonly,
            simplify=simplify,
            doit=doit,
            recursive=recursive,
            maxRecursionDepth=maxRecursionDepth
        )
    return res