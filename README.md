# critical-membrane-prewetting
  Simulations and Calculations for "Surface Densities Prewet a Near-Critical Membrane"
  	      https://doi.org/10.1101/2021.02.17.431700 
  Mason Rouches: mason.rouches@yale.edu
### Simulations
requires cython, numpy, networkx
#### Run simulations
	python run_sim_2phase.py # runs 2 phase simulation
	       Tether conc = 0.06
	       Bulk Chemical potential = -5.8
	       Ising temperature = 1.1Tc
	       Membrane Composition = 0.5
	python run_sim_3phase.py # runs 3 phase simulation
	       Tether conc = 0.08
	       Bulk Chemical potential = -4.5
	       Ising temperature = 0.85Tc
	       Membrane Composition = 0.35
### outputs 
        Creates directory: d[Dim1]-[Dim2]-J[bulk_coupling]_n[tether_conc]_l[poly_length]_u-[chem_potent]_hPhi[tether_coupling]_c[membrane_comp]_i[ising_coupling] where output files are stored in separate files for each coupling step
	   Each coupling step has 3 files - []_tethers.txt,[]_ising.txt,[]_polys.txt
	   Ising format is: iteration \t spins\n where spins are L^2 +/- 1 values for each position 
           Tether format is iteration \t tether_occ where tether_occ L^2 values 0/1 where 1 means tether is in position X
           Poly format is iteration \t sys \ t polytype \t poly_number \t positions \n
       	          where sys is 0 for system, 1 for reservoir
	    	  polytype is 0 for type 1, 1 for type 2
		  poly number counts the number of polymers in sys of poly type
		  positions is a poly_length*3 values where every 3 values corresponds to an (x,y,z) positon

### View Simulations
       generate_snapshots.py plots several random snapshots of each simulaion, along with averaged density profiles
       Subroutines for reading output files can be found There


### MFT Calculations
  Calculates the fixed-temperature phase diagram in terms of \lambda_{\rho},\lambda_{\psi} values. Translates to \psi,\rho values
  Run notebook sections sequentially
