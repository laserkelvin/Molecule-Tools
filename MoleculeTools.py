#!/bin/python


# MoleculeTools.py

from scipy import constants
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D

# Since I can't get PyBel to work properly, gonna try and
# write a few of my own routines to get the ball rolling

############################### Classes ###############################

# Molecule class, a dictionary of atoms. Might be a bit clunky to 
# reference...
class Molecule:
    Atoms = {}
    CalculatedCOM = False
    COM = np.array([0., 0., 0.])
    # Initialises by adding instances of atoms class to the dictionary
    # from XYZ format
    def __init__(self, InputXYZ):
        for atom in enumerate(InputXYZ):
            self.Atoms[str(atom[0])] = (Atom(atom[1][1], atom[1][2], atom[1][3], atom[1][0]))
    # Function to calculate the centre of mass for a given molecule in XYZ coordinates
    def CalculateCOM(self):
        CumSum = 0.                         # Cumulative sum of m_i * sum(r_i)
        TotalMass = 0.                      # Total mass of molecule
        for AtomNumber in range(len(self.Atoms)):               # Loop over all atoms
            CumSum = CumSum + self.Atoms[str(AtomNumber)].Mass * self.Atoms[str(AtomNumber)].Coordinates()
            TotalMass = TotalMass + self.Atoms[str(AtomNumber)].Mass
        self.CalculatedCOM = True           # Flag COM as having been calculated for plotting
        self.COM = (1. / TotalMass) * CumSum
        return self.COM
    # Function to plot up the molecule using an xyz matplotlib plot
    def Show(self):
        HydrogenRadius = 53.                          # in picometres
        AtomicRadii = {"H": 53. / HydrogenRadius,
                       "C": 67. / HydrogenRadius,
                       "O": 48. / HydrogenRadius}
        AtomicColours = {"H": "white",                # CPK colours
                         "C": "black",
                         "O": "red"}
        NAtoms = len(self.Atoms)
        X = np.zeros((NAtoms), dtype=float)                  # Arrays for holding xyz coordinates
        Y = np.zeros((NAtoms), dtype=float)
        Z = np.zeros((NAtoms), dtype=float)
        Colors = []
        Size = np.zeros((NAtoms), dtype=float)               # Size of the atoms
        for AtomNumber in range(len(self.Atoms)):               # Loop over all atoms
            X[AtomNumber] = self.Atoms[str(AtomNumber)].X
            Y[AtomNumber] = self.Atoms[str(AtomNumber)].Y
            Z[AtomNumber] = self.Atoms[str(AtomNumber)].Z
            AtomicSymbol = self.Atoms[str(AtomNumber)].Symbol
            Colors.append(AtomicColours[AtomicSymbol])          # work out the colour for atom
            Size[AtomNumber] = AtomicRadii[AtomicSymbol] * 150.
        fig = plt.figure()
        ax = plt.axes(projection = "3d")
        ax.scatter(X, Y, Z, s=Size, c=Colors)

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

############################### Moments of Inertia ###############################

# Here the routines take x,y,z as Angstroms, and mass in grams!
# XX diagonal elemnt of I
def IXX(y, z, mass):
    return mass * ((y * 1e-10)**2 + (z * 1e-10)**2)

# ZZ diagonal element of I
def IZZ(x, y, mass):
    return mass * ((x * 1e-10)**2 + (y * 1e-10)**2)

# YY diagonal element of I
def IYY(x, z, mass):
    return mass * ((x * 1e-10)**2 + (z * 1e-10)**2)

def IXY(x, y, mass):
    return mass * (x * y) * 1e-10

def IXZ(x, z, mass):
    return mass * (x * z) * 1e-10

def IYZ(y, z, mass):
    return mass * (y * z) * 1e-10

# Note on this: the standard orientation that Gaussian spits out is already rotated
# to the inertial frame; that means we only need to calculate the diagonal elements
# of the inertia tensor
def CalculateInertiaMatrix(Molecule):
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
    return constants.h / (8 * (np.pi**2) * (constants.c * 100) * Inertia)

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
