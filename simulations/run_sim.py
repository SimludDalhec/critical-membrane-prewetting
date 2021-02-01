#!/usr/bin/python
import sys
import polymer_membrane_sim


c_bound = float(sys.argv[1]) # Concentration of Tethers
Js = float(sys.argv[2])      # Bulk coupling to stop at
chem_potent = float(sys.argv[3]) # Chemical potential of bulk
hPhi = float(sys.argv[4])        # Bulk-Tether Coupling
comp = float(sys.argv[5])        # Membrane magnetization (1 = disordered,0 = ordered)
Ti = float(sys.argv[6])          # Ising coupling (Ti*Tc)

polymer_membrane_sim.main(c_bound,Js,chem_potent,hPhi,comp,Ti)
