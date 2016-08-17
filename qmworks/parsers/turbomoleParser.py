__author__ = "Felipe Zapata"

__all__ = ['readTurbomoleBasis', 'readTurbomoleMO']

# ===============> Standard libraries and third-party <=========================
import numpy as np
from pyparsing import (alphanums, alphas, delimitedList, Group, Literal,
                       nestedExpr, nums, OneOrMore, restOfLine, SkipTo,
                       Suppress, Word)

from qmworks.common import (AtomBasisData, AtomBasisKey, InfoMO)
from qmworks.parsers.parser import (floatNumber, floatNumberDot, natural)
from qmworks.utils import (concat, concatMap, zipWith, zipWith3)

from .cp2KParser import swapCoeff

# ===============> <===============
brackets = nestedExpr('[', ']', Word(alphanums))

dollar = Literal('$')

parenthesis = nestedExpr('(', ')', Word(alphanums))

pound = Suppress(Literal("#"))

slash = Literal('/')

openBra = Suppress(Literal('{'))

closeBra = Suppress(Literal('}'))

# ===========================> Basis File  <==================================

star = Suppress(Literal("*") + restOfLine)

header = dollar + restOfLine

parseKey = (Word(alphas, max=2)).setResultsName("atomLabel") + \
           (restOfLine).setResultsName("basisName")

basisFormat = delimitedList(Word(nums), '/')

contraction = Suppress(parenthesis + slash + brackets) + openBra + basisFormat + closeBra

basisHeader = natural + restOfLine

parseContr = pound + Suppress(Word(alphas, max=2)) + contraction

parseCoeff = Suppress(basisHeader) + OneOrMore(floatNumber)

parseBasisData = OneOrMore(Group(parseCoeff.setResultsName("contractions")))

parseBasis = (star + parseKey +
              parseContr.setResultsName("format") + star +
              parseBasisData.setResultsName("coeffs"))

topParseB = Suppress(header) + OneOrMore(Group(parseBasis))

# ==================> MOs <==================
headerMO = Suppress(SkipTo(Literal("[MO]")) + restOfLine)

sym = Literal("Sym") + restOfLine

spin = Literal("Spin") + restOfLine

occ = Literal("Occ") + restOfLine

numEntry = Suppress(natural + natural + Word(alphas) + natural + Word(alphanums)) + floatNumberDot

eigenValue = Suppress(sym + Literal("Ene=")) + floatNumberDot + Suppress(spin + occ)

eigenVector = OneOrMore(numEntry)

pair = eigenValue.setResultsName("eigenValue") + eigenVector.setResultsName("eigenVector")

topParserMO = Suppress(headerMO) + OneOrMore(Group(pair))
# =============================================================================
# Parsing From File


def readTurbomoleBasis(path):
    """Read Turbomole basis set"""
    bss = topParseB.parseFile(path)
    atoms = [xs.atomLabel.lower() for xs in bss]
    names = concat([xs.basisName.upper().split() for xs in bss])
    formats = [xs.format[:] for xs in bss]
    formats_int = map(lambda fss: [[int(x) for x in xs]
                                   for xs in fss], formats)
    rss = [rs.coeffs[:] for rs in bss]
    rawData = [[x.contractions[:] for x in rss[i]] for i in range(len(rss))]
    fst = lambda xs: xs[0]
    snd = lambda xs: xs[1]
    expos = list(map(mapFloat, [concatMap(fst, swapCoeff(2, rawData[i]))
                                for i in range(len(rawData))]))
    coeffs = list(map(mapFloat, [concatMap(snd, swapCoeff(2, rawData[i]))
                                 for i in range(len(rawData))]))
    basisData = zipWith(AtomBasisData)(expos)(coeffs)
    basiskey = zipWith3(AtomBasisKey)(atoms)(names)(formats_int)

    return basiskey, basisData


def readTurbomoleMO(path, nOrbitals, nOrbFuns=None):
    """
    read Turbomole molecular orbitals
    """
    rs = topParserMO.parseFile(path)
    es = np.array([safeFloat(r.eigenValue[0]) for r in rs])
    css = np.array([mapSafeFloat(r.eigenVector[:]) for r in rs])

    return InfoMO(es, css)

# Auxiliar functions


def mapFloat(xs):
    return list(map(float, xs))


def mapSafeFloat(xs):
    return list(map(safeFloat, xs))


def safeFloat(s):
    rs = s.split('.')
    if rs[0] == '-':
        return float('-0.' + rs[1])
    else:
        return float(s)
