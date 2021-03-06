CONFLIST HIL        HILBK HIL01 HIL02

NATOM    HILBK      6
NATOM    HIL01      11
NATOM    HIL02      11

IATOM    HILBK  N   0
IATOM    HILBK  H   1
IATOM    HILBK  CA  2
IATOM    HILBK  HA  3
IATOM    HILBK  C   4
IATOM    HILBK  O   5
IATOM    HIL01  CB  0
IATOM    HIL01 1HB  1
IATOM    HIL01 2HB  2
IATOM    HIL01  CG  3
IATOM    HIL01  ND1 4
IATOM    HIL01  CE1 5
IATOM    HIL01  HE1 6
IATOM    HIL01  NE2 7
IATOM    HIL01  HE2 8
IATOM    HIL01  CD2 9
IATOM    HIL01  HD2 10
IATOM    HIL02  CB  0
IATOM    HIL02 1HB  1
IATOM    HIL02 2HB  2
IATOM    HIL02  CG  3
IATOM    HIL02  ND1 4
IATOM    HIL02  HD1 5
IATOM    HIL02  CE1 6
IATOM    HIL02  HE1 7
IATOM    HIL02  NE2 8
IATOM    HIL02  CD2 9
IATOM    HIL02  HD2 10

ATOMNAME HILBK    0  N  
ATOMNAME HILBK    1  H  
ATOMNAME HILBK    2  CA 
ATOMNAME HILBK    3  HA 
ATOMNAME HILBK    4  C  
ATOMNAME HILBK    5  O  
ATOMNAME HIL01    0  CB 
ATOMNAME HIL01    1 1HB 
ATOMNAME HIL01    2 2HB 
ATOMNAME HIL01    3  CG 
ATOMNAME HIL01    4  ND1
ATOMNAME HIL01    5  CE1
ATOMNAME HIL01    6  HE1
ATOMNAME HIL01    7  NE2
ATOMNAME HIL01    8  HE2
ATOMNAME HIL01    9  CD2
ATOMNAME HIL01   10  HD2
ATOMNAME HIL02    0  CB 
ATOMNAME HIL02    1 1HB 
ATOMNAME HIL02    2 2HB 
ATOMNAME HIL02    3  CG 
ATOMNAME HIL02    4  ND1
ATOMNAME HIL02    5  HD1
ATOMNAME HIL02    6  CE1
ATOMNAME HIL02    7  HE1
ATOMNAME HIL02    8  NE2
ATOMNAME HIL02    9  CD2
ATOMNAME HIL02   10  HD2

#1.Basic Conformer Information: name, pka, em, rxn.
#Marilyn 6/11/03
#23456789A123456789B123456789C
PROTON   HIL01      0
PROTON   HIL02      0
PKA      HIL01      0.0
PKA      HIL02      0.0
ELECTRON HIL01      0
ELECTRON HIL02      0
EM       HIL01      0.0
EM       HIL02      0.0
RXN      HIL01      0.0
RXN      HIL02      0.0

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
CONNECT  HILBK  N   sp2       -1    C   0     CA  0     H
CONNECT  HILBK  H   s         0     N
CONNECT  HILBK  CA  sp3       0     N   0     C   0     CB  0     HA
CONNECT  HILBK  HA  s         0     CA
CONNECT  HILBK  C   sp2       0     CA  0     O   1     N
CONNECT  HILBK  O   sp2       0     C

CONNECT  HIL01  CB  sp3       0     CA  0     CG  0    1HB  0    2HB
CONNECT  HIL01 1HB  s         0     CB
CONNECT  HIL01 2HB  s         0     CB
CONNECT  HIL01  CG  sp2       0     CB  0     ND1 0     CD2
CONNECT  HIL01  ND1 sp2       0     CG  0     CE1
CONNECT  HIL01  CE1 sp2       0     ND1 0     NE2 0     HE1
CONNECT  HIL01  HE1 s         0     CE1
CONNECT  HIL01  NE2 sp2       0     CE1 0     CD2 0     HE2
CONNECT  HIL01  HE2 s         0     NE2
CONNECT  HIL01  CD2 sp2       0     NE2 0     CG  0     HD2
CONNECT  HIL01  HD2 s         0     CD2

CONNECT  HIL02  CB  sp3       0     CA  0     CG  0    1HB  0    2HB
CONNECT  HIL02 1HB  s         0     CB
CONNECT  HIL02 2HB  s         0     CB
CONNECT  HIL02  CG  sp2       0     CB  0     ND1 0     CD2
CONNECT  HIL02  ND1 sp2       0     CG  0     CE1 0     HD1
CONNECT  HIL02  HD1 s         0     ND1
CONNECT  HIL02  CE1 sp2       0     ND1 0     NE2 0     HE1
CONNECT  HIL02  HE1 s         0     CE1
CONNECT  HIL02  NE2 sp2       0     CE1 0     CD2 0
CONNECT  HIL02  CD2 sp2       0     NE2 0     CG  0     HD2
CONNECT  HIL02  HD2 s         0     CD2

#3.Atom Parameters: Partial Charges and Radii
# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   HIL    N   1.50
RADIUS   HIL    H   1.00
RADIUS   HIL    CA  2.00
RADIUS   HIL    HA  0.00
RADIUS   HIL    C   1.70
RADIUS   HIL    O   1.40
RADIUS   HIL    CB  2.00
RADIUS   HIL   1HB  0.00
RADIUS   HIL   2HB  0.00
RADIUS   HIL    CG  1.70
RADIUS   HIL    ND1 1.50
RADIUS   HIL    HD1 1.00
RADIUS   HIL    CE1 1.70
RADIUS   HIL    HE1 1.00
RADIUS   HIL    NE2 1.50
RADIUS   HIL    HE2 1.00
RADIUS   HIL    CD2 1.70
RADIUS   HIL    HD2 1.00

CHARGE   HILBK  N    -0.350
CHARGE   HILBK  H     0.250
CHARGE   HILBK  CA    0.100
CHARGE   HILBK  C     0.550
CHARGE   HILBK  O    -0.550

CHARGE   HIL01  CB    0.125
CHARGE   HIL01  CG    0.155
CHARGE   HIL01  NE2  -0.400
CHARGE   HIL01  HE2   0.400
CHARGE   HIL01  CE1   0.155
CHARGE   HIL01  HE1   0.125
CHARGE   HIL01  ND1  -0.560
CHARGE   HIL01  CD2  -0.125
CHARGE   HIL01  HD2   0.125

CHARGE   HIL02  CB    0.125
CHARGE   HIL02  CG    0.155
CHARGE   HIL02  NE2  -0.560
CHARGE   HIL02  CE1   0.155
CHARGE   HIL02  HE1   0.125
CHARGE   HIL02  ND1  -0.400
CHARGE   HIL02  HD1   0.400
CHARGE   HIL02  CD2  -0.125
CHARGE   HIL02  HD2   0.125
