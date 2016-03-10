#!/bin/python


# MoleculeTools.py

from scipy import constants
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from chemview import enable_notebook, MolecularViewer
enable_notebook()

# Since I can't get PyBel to work properly, gonna try and
# write a few of my own routines to get the ball rolling

############################### Classes ###############################

# Molecule class, a dictionary of atoms. Might be a bit clunky to 
# reference...
class Molecule:
    AtomScalingFactor = 400.
    # Initialises by adding instances of atoms class to the dictionary
    # from XYZ format
    def __init__(self, InputXYZ):
        self.COMCoordinates = False                             # Initial state is not in COM coordinates
        self.CalculatedCOM = False                              # Whether or not COM has been calculated already
        self.Atoms = {}                                         # Dictionary holding instances of Atom class
        self.COM = np.array([0., 0., 0.])                       # COM of molecule
        for atom in enumerate(InputXYZ):
            self.Atoms[str(atom[0])] = Atom(atom[1][1], atom[1][2], atom[1][3], atom[1][0])
        self.NAtoms = len(self.Atoms)
    # Function to calculate the centre of mass for a given molecule in XYZ coordinates
    def CalculateCOM(self):
        CumSum = 0.                                             # Cumulative sum of m_i * sum(r_i)
        TotalMass = 0.                                          # Total mass of molecule
        for AtomNumber in range(self.NAtoms):               # Loop over all atoms
            CumSum = CumSum + self.Atoms[str(AtomNumber)].Mass * self.Atoms[str(AtomNumber)].Coordinates()
            TotalMass = TotalMass + self.Atoms[str(AtomNumber)].Mass
        self.CalculatedCOM = True                               # Flag COM as having been calculated for plotting
        self.COM = (1. / TotalMass) * CumSum
        return self.COM
    # Function to shift the coordinates to centre of mass frame (I think it works)
    # by making the COM zero
    def Shift2COM(self):
        if self.COMCoordinates == True:
            print " Already in COM frame"
            pass
        else:
            self.CalculateCOM()
            COM = self.COM
            for AtomNumber in range(self.NAtoms):
                self.Atoms[str(AtomNumber)].X = self.Atoms[str(AtomNumber)].X + (-COM[0])
                self.Atoms[str(AtomNumber)].Y = self.Atoms[str(AtomNumber)].Y + (-COM[1])
                self.Atoms[str(AtomNumber)].Z = self.Atoms[str(AtomNumber)].Z + (-COM[2])
            self.COMCoordinates = True
    # Function that will format xyz
    def ExportXYZ(self):
        for AtomNumber in range(self.NAtoms):
            Coordinates = self.Atoms[str(AtomNumber)].Coordinates()
            print self.Atoms[str(AtomNumber)].Symbol + "\t" + str(Coordinates[0]) + "\t" + str(Coordinates[1]) + "\t" + str(Coordinates[2])

    # Function to plot up the molecule using an xyz matplotlib plot
    def Show(self):
        HydrogenRadius = 53.                                    # in picometres
        AtomicRadii = {"H": 53. / HydrogenRadius,
                       "C": 67. / HydrogenRadius,
                       "O": 48. / HydrogenRadius,
                       "COM": 50. / HydrogenRadius,
                       "X": 50. / HydrogenRadius}
        AtomicColours = {"H": "white",                          # CPK colours
                         "C": "black",
                         "O": "red",
                         "COM": "green",
                         "X": "purple"}
        NAtoms = len(self.Atoms)
        Colors = []
        if self.CalculatedCOM == False:
            X = np.zeros((NAtoms), dtype=float)                 # Arrays for holding xyz coordinates
            Y = np.zeros((NAtoms), dtype=float)
            Z = np.zeros((NAtoms), dtype=float)
            Size = np.zeros((NAtoms), dtype=float)              # Array for holding size of atoms
        else:
            X = np.zeros((NAtoms + 1), dtype=float)             # Arrays for holding xyz coordinates
            Y = np.zeros((NAtoms + 1), dtype=float)             # one more element for COM point
            Z = np.zeros((NAtoms + 1), dtype=float)
            Size = np.zeros((NAtoms + 1), dtype=float)          # Size of the atoms
        for AtomNumber in range(NAtoms):                        # Loop over all atoms
            X[AtomNumber] = self.Atoms[str(AtomNumber)].X
            Y[AtomNumber] = self.Atoms[str(AtomNumber)].Y
            Z[AtomNumber] = self.Atoms[str(AtomNumber)].Z
            AtomicSymbol = self.Atoms[str(AtomNumber)].Symbol
            Colors.append(AtomicColours[AtomicSymbol])          # work out the colour for atom
            Size[AtomNumber] = AtomicRadii[AtomicSymbol] * self.AtomScalingFactor
        if self.CalculatedCOM == True:                          # If we calculated COM before plot it too
            X[NAtoms] = self.COM[0]
            Y[NAtoms] = self.COM[1]
            Z[NAtoms] = self.COM[2]
            Size[NAtoms] = AtomicRadii["COM"] * self.AtomScalingFactor
            Colors.append(AtomicColours["COM"])
        fig = plt.figure()                                      # Use matplotlib to plot atoms in 3D scatter
        ax = plt.axes(projection = "3d")
        ax.w_xaxis.gridlines.set_lw(5.0)                        # This makes the grid lines thicker apparently
        ax.w_yaxis.gridlines.set_lw(5.0)
        ax.w_zaxis.gridlines.set_lw(5.0)
        ax.scatter(X, Y, Z, s=Size, c=Colors, depthshade=False)
        try:                                                    # if the principal axes have been calculated,
            type(self.PrincipalAxes)                            # show them on the plot too
            PX = np.zeros((3), dtype=float)                     # arrays for holding the principal axes
            PY = np.zeros((3), dtype=float)                     # components
            PZ = np.zeros((3), dtype=float)
            AxesColors = ["red", "green", "blue"]
            for AxisNumber in enumerate(self.PrincipalAxes):
                Vector = AxisNumber[1]
                ax.plot([0, Vector[0]], [0, Vector[1]], [0, Vector[2]],
                antialiased=True, color=AxesColors[AxisNumber[0]], linestyle="dashed")
        except AttributeError:
            pass

        plt.show()
    def GenerateBonds(self, Threshold=1.6):                     # Generate a list of bonds based on distance
        Bonds = []
        for AtomNumber1 in range(self.NAtoms):
            for AtomNumber2 in range(self.NAtoms):
                if AtomNumber1 == AtomNumber2:                  # If it's the same atom don't worry about it
                    pass
                else:
                    Distance = CalculateDistance(self.Atoms[str(AtomNumber1)], self.Atoms[str(AtomNumber2)])
                    if Distance <= Threshold:                   # If distance is less than threshold, it's a bond!
                        Bonds.append( (AtomNumber1, AtomNumber2) )
        self.Bonds = Bonds
        return Bonds
    def ChemView(self):                                         # Use ChemView module in notebook to plot molecule
        NAtoms = self.NAtoms
        Coordinates = np.zeros((NAtoms, 3), dtype=float)
        AtomSymbols = []
        for AtomNumber in range(NAtoms):                        # Make a list of XYZ and atomic symbol arrays
            Coordinates[AtomNumber] = self.Atoms[str(AtomNumber)].Coordinates()
            AtomSymbols.append(self.Atoms[str(AtomNumber)].Symbol)
        try:
            mv = MolecularViewer(Coordinates, topology={'atom_types': AtomSymbols,
                                                        'bonds': self.Bonds})
        except AttributeError:                                  # if self.Bonds isn't already generated, catch exception
            self.GenerateBonds()                                # If not already generated, make the bonds
            mv = MolecularViewer(Coordinates, topology={'atom_types': AtomSymbols,
                                                        'bonds': self.Bonds})
        mv.ball_and_sticks()
        return mv
    # Routine effectively copied from down below, but adapted for use as molecule method
    # Calculates the moment of inertia matrix without assuming that it's diagonal...
    def CalculateInertiaMatrix(self):
        """
        Calculates the Inertia Tensor assuming a center-of-mass frame.
        """
        if self.COMCoordinates == False:                        # come up with warning about not being in COM frame
            print " Warning: not in COM frame! Inertia matrix will be wrong!"
        else:
            pass                                                # no need to worry if already in COM frame
        InertiaMatrix = np.zeros((3,3), dtype=float)
        self.InertiaMatrix = InertiaMatrix                      # Zero the matrix beforehand
        for AtomNumber in range(self.NAtoms):                   # Loop over atoms in molecule
            Coordinates = self.Atoms[str(AtomNumber)].Coordinates() * 1e-10      # Retrieve coordinates, convert to metres
            Mass = self.Atoms[str(AtomNumber)].Mass                              # Retrieve mass of atom
            # Work out diagonal elements of the matrix
            InertiaMatrix[0,0] = InertiaMatrix[0,0] + IXX(Coordinates[1], Coordinates[2], Mass)
            InertiaMatrix[1,1] = InertiaMatrix[1,1] + IYY(Coordinates[0], Coordinates[2], Mass)
            InertiaMatrix[2,2] = InertiaMatrix[2,2] + IZZ(Coordinates[0], Coordinates[1], Mass)
            # Calculate off-diagonal elements of matrix
            InertiaMatrix[0,1] = InertiaMatrix[0,1] + IXY(Coordinates[0], Coordinates[1], Mass)
            InertiaMatrix[0,2] = InertiaMatrix[0,2] + IXZ(Coordinates[0], Coordinates[2], Mass)
            InertiaMatrix[1,2] = InertiaMatrix[1,2] + IYZ(Coordinates[1], Coordinates[2], Mass)
        # Apply sign change
        InertiaMatrix[0,1] = -InertiaMatrix[0,1]
        InertiaMatrix[0,2] = -InertiaMatrix[0,2]
        InertiaMatrix[1,2] = -InertiaMatrix[1,2]
        # Symmetrize 
        self.InertiaMatrix = (InertiaMatrix + InertiaMatrix.T) / 2.
        return self.InertiaMatrix
    # Diagonalise the moments of inerta matrix and work out the principal moments of inertia
    def PrincipalMoments(self):
        if type(self.InertiaMatrix) == None:                    # check if MOI matrix has been calculated
            self.CalculateInertiaMatrix()
        else:
            Diagonal = np.linalg.eig(self.InertiaMatrix)     # Ignore the eigenvectors
            Eigenvalues = Diagonal[0]
            Eigenvectors = Diagonal[1]
            self.PMI = PMI2ABC(Eigenvalues)                  # Return the rotational constants in 1/cm
            self.PrincipalAxes = Eigenvectors                # Return the eigenvectors for plotting
            return self.PMI
    def SumMass(self):                                          
        Mass = 0.
        for AtomNumber in range(self.NAtoms):
            Mass = Mass + self.Atoms[str(AtomNumber)].Mass
        self.Mass = Mass
        return self.Mass

# Atom class, has attributes of the xyz coordinates as well as its symbol and mass
class Atom:
    X = 0.
    Y = 0.
    Z = 0.
    Symbol = " "
    Mass = 0.
    def __init__(self, X, Y, Z, Symbol):
        self.X = float(X)
        self.Y = float(Y)
        self.Z = float(Z)
        self.Symbol = Symbol
        self.Mass = Symbol2Mass(self.Symbol)
    # Returns the coordinates of an atom
    def Coordinates(self):
        return np.array([self.X, self.Y, self.Z])

############################### I/O ###############################

# Function to read the output xyz coordinates. The below functions
# will work better if this is rotated to the inertial frame of reference!
# i.e. standard orientation in Gaussian
def ReadXYZ(File):
    f = open(File, "r")
    fc = f.readlines()
    f.close()
    NAtoms = len(fc)
    Coordinates = []
    for line in range(NAtoms):
        Coordinates.append(fc[line].split())
    return Coordinates

def ReadNormalModes(File):
    f = open(File, "r")
    fc = f.readlines()
    f.close()
    NModes = len(fc)
    NormalModes = []
    for line in range(NModes):
        NormalModes.append(fc[line].split())
    return NormalModes

############################### Tools ###############################

# Function that stores a library with all the atomic masses
# and returns it for whatever atom you specify
def Symbol2Mass(Atom):
    # Masses are taken from:
    # http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some
    MassList = {"H": 1.00782503223,
                "C": 12.00,
                "O": 15.99491461957,
                "X": 0.0                                    # Dummy atom
                }
#    return MassList[Atom]                                  # Returns mass in amu
    return (MassList[Atom] / constants.Avogadro) / 1e3      # Return mass in kg

# Calculates the distance between two atoms
def CalculateDistance(A, B):
    dX = (A.X - B.X)**2
    dY = (A.Y - B.Y)**2
    dZ = (A.Z - B.Z)**2
    return np.sqrt(dX + dY + dZ)

# Uses cartesian points A, B, C to calculate the angle
# formed by A - B - C using vectors. Atom B is the centre of the angle!
def CalculateAngle(A, B, C):
    AB = B.Coordinates() - A.Coordinates()           # Calculate vector AB
    BC = C.Coordinates() - B.Coordinates()           # Calculate vector BC
    ABLength = CalculateDistance(A, B)               # Work out magnitude of AB
    BCLength = CalculateDistance(B, C)               # Magnitude of BC
    DotProduct = np.dot(AB, BC)                      # Dot product of AB dot BC
    # Return the angle formed by A - B - C in degrees
    return 180. - (np.arccos(DotProduct / (ABLength * BCLength)) * (180. / np.pi))

# Function to return the reduced mass of fragments
# Takes a list of masses in whatever units
def CalculateReducedMass(Masses):
    ReducedMass = 0.
    for mass in enumerate(Masses):
        ReducedMass = ReducedMass + (1. / mass[1])
    return 1. / ReducedMass

############################### Moments of Inertia ###############################

# Here the routines take x,y,z as Angstroms, and mass in kilograms!
# XX diagonal elemnt of I
def IXX(y, z, mass):
    return mass * ((y)**2 + (z)**2)

# ZZ diagonal element of I
def IZZ(x, y, mass):
    return mass * ((x)**2 + (y)**2)

# YY diagonal element of I
def IYY(x, z, mass):
    return mass * ((x)**2 + (z)**2)
def IXY(x, y, mass):
    return mass * (x * y)

def IXZ(x, z, mass):
    return mass * (x * z)

def IYZ(y, z, mass):
    return mass * (y * z)

# Note on this: the standard orientation that Gaussian spits out is already rotated
# to the inertial frame; that means we only need to calculate the diagonal elements
# of the inertia tensor
def OldCalculateInertiaMatrix(Molecule):
    InertiaMatrix = np.zeros((3,3),dtype=float)
    for atom in enumerate(Molecule):
        # Calculate diagonal elements of inertia matrix
        # Enumerate loops through each atom of the molecule;
        # indexes 0:mass, 1:x, 2:y, 3:z
        InertiaMatrix[0,0] = InertiaMatrix[0,0] + IXX(float(atom[1][2]), float(atom[1][3]), Symbol2Mass(atom[1][0]))
        InertiaMatrix[1,1] = InertiaMatrix[1,1] + IYY(float(atom[1][1]), float(atom[1][3]), Symbol2Mass(atom[1][0]))
        InertiaMatrix[2,2] = InertiaMatrix[2,2] + IZZ(float(atom[1][1]), float(atom[1][2]), Symbol2Mass(atom[1][0]))
        # Calculate triangle of off-diagonal elements of inertia matrix
#        InertiaMatrix[0,1] = InertiaMatrix[0,1] + IXY(float(atom[1][1]), float(atom[1][2]), Symbol2Mass(atom[1][0]))
#        InertiaMatrix[0,2] = InertiaMatrix[0,2] + IXZ(float(atom[1][1]), float(atom[1][3]), Symbol2Mass(atom[1][0]))
#        InertiaMatrix[1,2] = InertiaMatrix[1,2] + IYZ(float(atom[1][2]), float(atom[1][3]), Symbol2Mass(atom[1][0]))
        # Symmetrise the matrix
    return (InertiaMatrix + InertiaMatrix.T) / 2.

# Converts the principle moments into rotational constants in 1/cm
# These should agree with the rotational constants in Gaussian
def PMI2ABC(Inertia):
    return constants.h / (8 * (np.pi)**2 * (constants.c * 100.) * Inertia)

############################### Normal mode analysis ###############################

# Take an input cartesian coordinate (x,y,z) as well as the corresponding
# normal mode displacement for that x,y,z
def CalculateDisplacement(x, a):
    return (a * x)**2

# Function to calculate displaced geometries using normal mode displacements
def NormalMode2Cartesian(Molecule, NormalModes, NormalDisplacment):
    DisplacedGeometry = np.zeros((len(Molecule), 3), dtype=float)
    for atom in enumerate(Molecule):
        Displacement = NormalModes[atom[0]] 
