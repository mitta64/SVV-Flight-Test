# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 09:28:53 2020

@author: Thomas
"""

#Imports
import numpy as np
from Cit_par import *
from AnalyticalModel import ShortPeriodSimplified, PhugoidSimplified, HeavilyDampedApriodicRollSimplified, DutchRollSimplified, ApriodicSpiralSimplified, DutchandAperiodRollSimplified


#Functions to calculate certain eigenvalue parameters
def Thalf(eigreal):
    Thalf = (np.log(0.5)/eigreal)*(b/V0)
    return Thalf

def Period(eigimag):
    P = (2*np.pi/eigimag)*(b/V0)
    return P

def Dampratio(eigreal, eigimag):
    damp = -eigreal / (np.sqrt(eigreal**2 + eigimag**2))
    return damp

#The eigenvalues determined by the numerical model
eigsym = np.array([-1.17057194+1.34384154j, -1.17057194-1.34384154j, -0.00782837+0.14531042j, -0.00782837-0.14531042j])
eigasym = np.array([-3.793, -0.2556 +2.009j, -0.2556 -2.009j, 0.01054])


###symmetric eigenmotions###
#Short Period
Lsp1 = eigsym[0]
Lsp2 = eigsym[1]

Thalfsp = Thalf(Lsp1.real)
Psp = Period(Lsp1.imag)
zetasp = Dampratio(Lsp1.real, Lsp1.imag)
#print(Lsp1, Lsp2, Thalfsp, Psp, zetasp)
#print(ShortPeriodSimplified())

    #The dampratio seems to be exactly the same, makes me think that the magintudes are different but ratios betwee the real and imag part are almost identical
Lsp1ver = ShortPeriodSimplified()[0]
Lsp2ver = ShortPeriodSimplified()[1]
#print(Lsp1/Lsp1.real, Lsp1ver/Lsp1ver.real )
#print(Lsp1/Lsp1.real, Lsp2/Lsp2.real, Thalf((Lsp1/Lsp1.real).real), Period((Lsp1/Lsp1.real).imag), Dampratio((Lsp1/Lsp1.real).real, (Lsp1/Lsp1.real).imag))
#print(Lsp1ver/Lsp1ver.real, Lsp2ver/Lsp1ver.real, Thalf((Lsp2ver/Lsp1ver.real).real) , Period((Lsp2ver/Lsp1ver.real).imag), Dampratio((Lsp2ver/Lsp1ver.real).real, (Lsp2ver/Lsp1ver.real.imag)))

#Phugoid
Lph1 = eigsym[2]
Lph2 = eigsym[3]

Thalfph = Thalf(Lph1.real)
Pph = Period(Lph1.imag)
zetaph = Dampratio(Lph1.real, Lph1.imag)

#print(Lph1, Lph2, Thalfph, Pph, zetaph)
#print(PhugoidSimplified())
#print()

#looking at the damping error
error = ((0.0441 - 0.1034)/0.0441)*100
zetaana = PhugoidSimplified()[-1]
errorph = ((zetaph - zetaana)/zetaana)*100
#print(error, errorph)

#Ratio real and imag parts
ratioph = Lph2.real/Lph2.imag
ratiover = (PhugoidSimplified()[1]).real/(PhugoidSimplified()[1]).imag
#print(ratioph, ratiover)
#The dampratio seems to be exactly the same, makes me think that the magintudes are different but ratios betwee the real and imag part are almost identical
#Lph1ver = PhugoidSimplified()[0]
#Lph2ver = PhugoidSimplified()[1]
#print(Lph1/Lph1.real, Lph1ver/Lph1ver.real)
#print(Lph1/Lph1.imag, Lph2/Lph2.imag, Thalf((Lph1/Lph1.imag).real), Period((Lph1/Lph1.imag).imag), Dampratio((Lph1/Lph1.imag).real, (Lph1/Lph1.imag).imag))
#print(Lph1ver/Lph1.imag, Lph2ver/Lph1.imag, Thalf((Lph2ver/Lph1.imag).real) , Period((Lph2ver/Lph1.imag).imag), Dampratio((Lph2ver/Lph1.imag).real, (Lph2ver/Lph1.imag).imag))



###Asymmetric eigenmotions###
#Heavily damped aperiodic motion
Lapr = eigasym[0]
Thalfaperroll = Thalf(Lapr)
print(Lapr, Thalfaperroll)
print(HeavilyDampedApriodicRollSimplified())
print('error', ((-0.3291 + 0.4629)/-0.3291)*100, ((Lapr + HeavilyDampedApriodicRollSimplified()[0])/Lapr)*100)
print()

#Dutch roll
Ldr1 = eigasym[1]
Ldr2 = eigasym[2]
Thalfdr = Thalf(Ldr1.real)
Pdr = Period(Ldr1.imag)
zetadr = Dampratio(Ldr1.real, Ldr1.imag)
print(Ldr1, Ldr2, Thalfdr, Pdr, zetadr)
print(DutchandAperiodRollSimplified()[1], DutchandAperiodRollSimplified()[2], DutchandAperiodRollSimplified()[4], DutchandAperiodRollSimplified()[5], DutchandAperiodRollSimplified()[-1])
print()

#Aperidic spiral
Lspir = eigasym[-1]
Thalfspir = Thalf(Lspir)
print(Lspir, Thalfspir)
print(ApriodicSpiralSimplified())







