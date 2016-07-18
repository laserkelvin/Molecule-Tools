# CoulombMatrix.py


""" Routines for calculating the Coulomb Matrix
    representation of a molecule.

"""

import numpy as np
import glob
import pandas as pd

#############################################################################

""" Coulomb Matrix routines """

#############################################################################

def CMDiagonal(Charge):
    """ Diagonal elements only depend on charge of one atom """
    return 0.5 * Charge**(2.4)

def CMOffDiagonal(Atom1, Atom2):
    """ Off-diagonal elements of Coulomb Matrix from Hansen et al. JCTC, 2013 """
    Denominator = np.sqrt(np.sum((Atom1.GetCoordinates() - Atom2.GetCoordinates())**2.))
    Numerator = Atom1.GetCharge() * Atom2.GetCharge()
    return Numerator / Denominator

def CalculateCoulombMatrix(Molecule):
    """ For a given dictionary input of atom classes,
        return the associated Coulomb matrix
    """
    NAtoms = len(Molecule)          # Number of atoms in the molecule
    CoulombMatrix = np.zeros((NAtoms, NAtoms), dtype=float)   # initialise array
    for i in xrange(NAtoms):        # Loop over atoms i
        for j in xrange(NAtoms):    # Loop over atoms j
            if i == j:              # Calculate diagonal element as effective charge
                CoulombMatrix[i, j] = CMDiagonal(Molecule[i].GetCharge())
            else:                   # Calculate off-diagonal elements
                CoulombMatrix[i, j] = CMOffDiagonal(Molecule[i], Molecule[j])
    return CoulombMatrix

#############################################################################

""" Parsing, molecule and atom routines """

#############################################################################

class __Atom__:
    def __init__(self, Coordinate, Symbol):
        """ Label is atom symbol string,
            Coordinates is 3-tuple list
            
            Charges are taken from this list:
            https://en.wikipedia.org/wiki/Effective_nuclear_charge#Values
        """
        ChargeTable = {"C": 12.026,
                       "B": 9.677,
                       "H": 1.,
                       "O": 16.603,
                       "N": 14.346,
                       "He": 1.688,
                       "X": 0.0,               # Dummy atom
                      }
        self.Coordinates = np.array(Coordinate)
        self.Symbol = Symbol
        self.Charge = ChargeTable[self.Symbol]
    # Returns the coordinates of an atom
    def GetCoordinates(self):
        return self.Coordinates
    
    def GetSymbol(self):
        return self.Symbol
    
    def GetCharge(self):
        return self.Charge

def MoleculeFromXYZ(File):
    """ Read in a chemical file format .xyz """
    f = open(File, "r")
    fc = f.readlines()[2:]            # Skip the number of atoms and comment line
    f.close()
    NAtoms = len(fc)
    Molecule = dict()
    for Line in range(NAtoms):
        SplitLine = fc[Line].split()
        Symbol = SplitLine[0]                                 # First item is atom symbol
        Coordinates = np.array([SplitLine[1], SplitLine[2], SplitLine[3]])
        Coordinates = Coordinates.astype(np.float)            # convert coordinates to numpy float array
        Molecule[Line] = __Atom__(Coordinates, Symbol)            # Populate dictionary
    return Molecule

def StripExtension(Filename):
    """ Strips the extension of the file only. """
    return Filename.split(".")[0]

def StripFolder(Filename):
    """ Strips the path information from file """
    return Filename.split("/")[-1]

#############################################################################

""" Pre-processing routines """

#############################################################################

def PadDummyAtoms(Molecule, LargestMolecule):
    """ Will pad a Molecule with dummy atoms based on
        the largest number of atoms
    """
    for Atom in range(len(Molecule), LargestMolecule):
            Molecule[Atom] = __Atom__(np.random.rand(3), "X")  # random numbers to prevent division by zero
    return Molecule

def BatchCalculate(Folder=None):
    """ Takes a string input for the folder, and will
        look in the relative path for all files containing
        .xyz extension. 

        If nothing is specified, it'll look in the current folder.

        Parses every .xyz file in folder, reads in the atoms
        and calculates the eigenvalues of the Coulomb matrix.

        The return data is a DataFrame with keys as the file name
        and values as the principal moments of the Coulomb matrix.

        This way we solve the two problems outlined in:
        K. Hansen et al. JCTC 130730111426003 (2013).

        1. Different dimensionality of molecules:
            Solved by padding molecules with dummy atoms

        2. Descriptors of molecules:
            We're using the eigenvalues of the Coulomb matrix.
            Should be quite different for two similar molecules,
            but I need to test this more.
    """
    if Folder == None:
        Path = "./*.xyz"
    else:
        Path = Folder + "/*.xyz"

    """ Initialise some data """
    FileList = glob.glob(Path)
    LargestMolecule = 0           # Stores the number of atoms of the largest molecule
    DataFrame = pd.DataFrame()    # Dataframe that holds the eigenvalues

    """ First pass figures out how many atoms the largest molecule has """
    for File in FileList:
        Molecule = MoleculeFromXYZ(File)
        if len(Molecule) > LargestMolecule:               # if this is larger than the current largest
            LargestMolecule = len(Molecule)               # set it as the largest number of atoms
            LargestName = StripExtension(File)
    print "The largest molecule has " + str(LargestMolecule) + " atoms."
    print "The molecule is " + LargestName

    """ Second pass actually calculates the Coulomb matrices """
    for File in FileList:
        Molecule = MoleculeFromXYZ(File)
        if len(Molecule) < LargestMolecule:
            Molecule = PadDummyAtoms(Molecule, LargestMolecule)
        else:
            pass
        Name = StripFolder(File)                          # Strip the extension and path
        Name = StripExtension(Name)
        CoulombMatrix = CalculateCoulombMatrix(Molecule)  # Calculate the Coulomb matrix
        Eigenvalues = np.linalg.eigvals(CoulombMatrix)    # Calculate eigenvalues of CM
        DataFrame[Name] = Eigenvalues
    return DataFrame