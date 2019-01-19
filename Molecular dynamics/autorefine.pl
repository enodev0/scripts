#!/usr/bin/env perl
use strict;
use warnings;

# This script tries to automate a basic protein-in-water MD simulation in GROMACS, result analysis
# included. This is mostly useful for automated tertiary structure refinement purposes.
# Comment out the halt_and_prompt() function call from the main loop below if interactive
# use is not desired.

# The steps followed here are based on: http://www.mdtutorials.com/gmx/lysozyme/index.html

# NOTE: This script assumes the PDB file is of sound quality, and that you want to use the
# OPLS-AA all-atom forcefield to model your protein.
# Certain parameters, particularly the non-bonded interaction cutoffs (eg: van der waals cutoff)
# need to be modified according to your forcefield of choice.

# Requirement: GROMACS v2016.1


# 1ps = 1000fs
# 1ns = 1000ps = 10^6fs

################################################################################
# 	                        PARAMETERS                                     #
################################################################################
my $PROTEIN_NAME 	      =   "lcfa_synthetase";
my $PROTEIN_ID                =   "facs";
my $SOLVENT_MODEL             =   "spc216"; # generic 3-point water model
my $SOLVENT_TYPE              =   "spce";   # water
my $BOX_MOLECULE_CONTINGENCY  =   0.75;      # nm
my $BOX_TYPE                  =   "cubic";  # dodecahedron?
my $ION_CONCENTRATION         =   0.04;      # Molar of NaCl, after automatic charge neutralization

my $STEPWISE_HALT             =   0;        # false
my $INTEGRATOR_STEEP          =   "steep";
my $INTEGRATOR_LEAPFROG       =   "md";
my $UNIT_CLOCKTICK            =   0.002;    # 2 femtoseconds, written as picoseconds
my $MDRUN_DURATION_NANO       =   10;        # in nanoseconds, simulation duration
my $NPT_DURATION_PICO         =   30;      # in picoseconds, Number-Pressure-Temperature
my $NVT_DURATION_PICO         =   10;      # in picoseconds, Number-Volume-Temperature
my $ION_GENERATION_PICO       =   50;      # in picoseconds
my $ENERGY_MINIM_PICO         =   50;      # in picoseconds
my $NPT_STEPCOUNT             =   ($NPT_DURATION_PICO * 1000)/($UNIT_CLOCKTICK * 1000);
my $NVT_STEPCOUNT             =   ($NVT_DURATION_PICO * 1000)/($UNIT_CLOCKTICK * 1000);
my $MDRUN_STEPCOUNT           =   ($MDRUN_DURATION_NANO * 1000000)/($UNIT_CLOCKTICK * 1000);
my $ION_GEN_STEPCOUNT         =   ($ION_GENERATION_PICO * 100);
my $ENR_MIN_STEPCOUNT         =   ($ENERGY_MINIM_PICO * 100);
my $UNIT_CLOCKTICK_FEMTO      =   ($UNIT_CLOCKTICK * 1000);

################################################################################
#                              COMMAND LIST                                    #
################################################################################
my $CMD_PDB2GMX      = "echo '10' | gmx pdb2gmx -f $PROTEIN_ID.pdb -o conf.gro -water $SOLVENT_TYPE"; # select opls-aa forcefield
my $CMD_EDITCONF     = "gmx editconf -f conf.gro -o box.gro -c -d $BOX_MOLECULE_CONTINGENCY -bt $BOX_TYPE";
my $CMD_SOLVATE      = "gmx solvate -cp box.gro -cs $SOLVENT_MODEL.gro -o solvated.gro -p topol.top";
my $CMD_GROMPION     = "gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr";
my $CMD_GENION       = "echo '13' | gmx genion -s ions.tpr -neutral -conc $ION_CONCENTRATION -p topol.top -o ions.gro"; # choose 13 SOL
my $CMD_GROMPPOT     = "gmx grompp -f em.mdp -c solvated.gro -p topol.top -o em.tpr";
my $CMD_PE_MINIMIZE  = "gmx mdrun -v -deffnm em";
my $CMD_PLOT_PE      = "echo '10 0' | gmx energy -f em.edr -o potential.xvg"; # select 10 0
my $CMD_GROMPNVT     = "gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr";
my $CMD_NVT_EQUIB    = "gmx mdrun -deffnm nvt";
my $CMD_PLOT_NVT     = "echo '15 0' | gmx energy -f nvt.edr -o temperature.xvg"; # select 15 0
my $CMD_GROMPNPT     = "gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr";
my $CMD_NPT_EQUIB    = "gmx mdrun -deffnm npt";
my $CMD_PLOT_NPT     = "echo '16 0' | gmx energy -f npt.edr -o pressure.xvg"; # select 16 0
my $CMD_PLOT_DEN     = "echo '22 0' | gmx energy -f npt.edr -o density.xvg";  # select 22 0
my $CMD_GROMPMD      = "gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o simulation.tpr";
my $CMD_SIMULATE     = "gmx mdrun -deffnm simulation";
my $CMD_TRJCONV      = "echo '0' | gmx trjconv -s simulation.tpr -f simulation.xtc -o simulation-noPBC.xtc -pbc mol -ur compact"; # choose System
my $CMD_PLOT_RMSD    = "echo '4' | gmx rms -s simulation.tpr -f simulation-noPBC.xtc -o rmsd.xvg -tu ns"; # Choose Backbone
my $CMD_PLOT_GYRA    = "echo '4' | gmx gyrate -s simulation.tpr -f simulation-noPBC.xtc -o gyrate.xvg";  # Choose Backbone
my $CMD_PLOT_RMSF    = "echo '3' | gmx rmsf -s simulation.tpr -f simulation.xtc -o rmsf.xvg -oq bfac.pdb"; # select C-alpha (3)
################################################################################
#                                 MDRUN FILES                                  #
################################################################################
# em.mdp
my $EMI_MDP = "; em.mdp - used as input into grompp to generate em.tpr
integrator	    = $INTEGRATOR_STEEP      ; Algorithm (steep = steepest descent minimization)
emtol		    = 1000.0  	             ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      	    = 0.01                   ; Energy step size
nsteps		    = $ENR_MIN_STEPCOUNT     ; Maximum number of (minimization) steps to perform (default=50000)

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		    	     ; Frequency to update the neighbor list and long range forces
cutoff-scheme  	    = Verlet                 ; REMOVE THIS LINE IF USING GPU
ns_type		    = grid		     ; Method to determine neighbor list (simple, grid)
coulombtype	    = PME		     ; Treatment of long range electrostatic interactions
rcoulomb	    = 1.0		     ; Short-range electrostatic cut-off
rvdw		    = 1.0		     ; Short-range Van der Waals cut-off
pbc		    = xyz 		     ; Periodic Boundary Conditions (yes/no)";
############################################################
# ions.mdp
############################################################
my $ION_MDP = "; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator	    = $INTEGRATOR_STEEP	     ; Algorithm (steep = steepest descent minimization)
emtol		    = 1000.0  	             ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep              = 0.01      	     ; Energy step size
nsteps		    = $ION_GEN_STEPCOUNT     ; Maximum number of (minimization) steps to perform (default=50000)

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		    	     ; Frequency to update the neighbor list and long range forces
cutoff-scheme       = Verlet                 ; REMOVE THIS LINE IF USING GPU
ns_type		    = grid		     ; Method to determine neighbor list (simple, grid)
coulombtype	    = PME		     ; Treatment of long range electrostatic interactions
rcoulomb	    = 1.0		     ; Short-range electrostatic cut-off
rvdw		    = 1.0		     ; Short-range Van der Waals cut-off
pbc		    = xyz 		     ; Periodic Boundary Conditions (yes/no)";
############################################################
# npt.mdp
############################################################
my $NPT_MDP = "title	= $PROTEIN_NAME NPT equilibration
define		        = -DPOSRES	         ; position restrain the protein
; Run parameters
integrator	        = $INTEGRATOR_LEAPFROG	 ; leap-frog integrator
nsteps		        = $NPT_STEPCOUNT	 ; ($UNIT_CLOCKTICK_FEMTO * $NPT_STEPCOUNT) fs = $NPT_DURATION_PICO ps
dt		        = $UNIT_CLOCKTICK	 ; $UNIT_CLOCKTICK_FEMTO fs
; Output control
nstxout		        = 500		         ; save coordinates every 1.0 ps
nstvout		        = 500		         ; save velocities every 1.0 ps
nstenergy	        = 500		         ; save energies every 1.0 ps
nstlog		        = 500		         ; update log file every 1.0 ps
; Bond parameters
continuation	        = yes		         ; Restarting after NVT
constraint_algorithm    = lincs	    	         ; holonomic constraints
constraints	        = all-bonds	         ; all bonds (even heavy atom-H bonds) constrained
lincs_iter	        = 1	                 ; accuracy of LINCS
lincs_order	        = 4		         ; also related to accuracy
; Neighborsearching
cutoff-scheme           = Verlet                 ; REMOVE THIS LINE IF USING GPU
ns_type		        = grid		         ; search neighboring grid cells
nstlist		        = 10	    	         ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb	        = 1.0		         ; short-range electrostatic cutoff (in nm)
rvdw		        = 1.0		         ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype	        = PME		         ; Particle Mesh Ewald for long-range electrostatics
pme_order	        = 4		         ; cubic interpolation
fourierspacing	        = 0.16		         ; grid spacing for FFT
; Temperature coupling is on
tcoupl		        = V-rescale	         ; modified Berendsen thermostat
tc-grps		        = Protein Non-Protein    ; two coupling groups - more accurate
tau_t		        = 0.1	  0.1	         ; time constant, in ps
ref_t		        = 300 	  300	         ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl		        = Parrinello-Rahman      ; Pressure coupling on in NPT
pcoupltype	        = isotropic	         ; uniform scaling of box vectors
tau_p		        = 2.0		         ; time constant, in ps
ref_p		        = 1.0		         ; reference pressure, in bar
compressibility         = 4.5e-5	         ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com
; Periodic boundary conditions
pbc		        = xyz		         ; 3-D PBC
; Dispersion correction
DispCorr	        = EnerPres	         ; account for cut-off vdW scheme
; Velocity generation
gen_vel		        = no		         ; Velocity generation is off ";
############################################################
# nvt.mdp
############################################################
my $NVT_MDP = "title    = $PROTEIN_NAME NVT equilibration
define		        = -DPOSRES	         ; position restrain the protein
; Run parameters
integrator	        = $INTEGRATOR_LEAPFROG	 ; leap-frog integrator
nsteps		        = $NVT_STEPCOUNT	 ; ($UNIT_CLOCKTICK_FEMTO * $NVT_STEPCOUNT) fs = $NVT_DURATION_PICO ps
dt		        = $UNIT_CLOCKTICK	 ; $UNIT_CLOCKTICK_FEMTO fs
; Output control
nstxout		        = 500		         ; save coordinates every 1.0 ps
nstvout		        = 500		         ; save velocities every 1.0 ps
nstenergy	        = 500		         ; save energies every 1.0 ps
nstlog		        = 500		         ; update log file every 1.0 ps
; Bond parameters
continuation	        = no		         ; first dynamics run
constraint_algorithm    = lincs	                 ; holonomic constraints
constraints	        = all-bonds	         ; all bonds (even heavy atom-H bonds) constrained
lincs_iter	        = 1		         ; accuracy of LINCS
lincs_order	        = 4		         ; also related to accuracy
; Neighborsearching
cutoff-scheme           = Verlet                 ; REMOVE THIS LINE IF USING GPU
ns_type		        = grid		         ; search neighboring grid cells
nstlist		        = 10		         ; 20 fs, largely irrelevant with Verlet
rcoulomb	        = 1.0		         ; short-range electrostatic cutoff (in nm)
rvdw		        = 1.0		         ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype	        = PME	                 ; Particle Mesh Ewald for long-range electrostatics
pme_order	        = 4		         ; cubic interpolation
fourierspacing	        = 0.16		         ; grid spacing for FFT
; Temperature coupling is on
tcoupl		        = V-rescale	         ; modified Berendsen thermostat
tc-grps		        = Protein Non-Protein	 ; two coupling groups - more accurate
tau_t			= 0.1	  0.1            ; time constant, in ps
ref_t			= 300 	  300            ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl		        = no 		         ; no pressure coupling in NVT
; Periodic boundary conditions
pbc		        = xyz		         ; 3-D PBC
; Dispersion correction
DispCorr	        = EnerPres	         ; account for cut-off vdW scheme
; Velocity generation
gen_vel		        = yes		         ; assign velocities from Maxwell distribution
gen_temp	        = 300		         ; temperature for Maxwell distribution
gen_seed	        = -1		         ; generate a random seed";
############################################################
# md.mdp
############################################################
my $MDR_MDP = "title	= $PROTEIN_NAME MD simulation
; Run parameters
integrator	        = $INTEGRATOR_LEAPFROG	 ; leap-frog integrator
nsteps		        = $MDRUN_STEPCOUNT	 ; ($UNIT_CLOCKTICK_FEMTO * $MDRUN_STEPCOUNT) fs = $MDRUN_DURATION_NANO ns
dt		        = $UNIT_CLOCKTICK	 ; $UNIT_CLOCKTICK_FEMTO fs
; Output control
nstxout		        = 5000		         ; save coordinates every 10.0 ps
nstvout		        = 5000		         ; save velocities every 10.0 ps
nstenergy	        = 5000		         ; save energies every 10.0 ps
nstlog		        = 5000		         ; update log file every 10.0 ps
nstxout-compressed      = 5000                   ; save compressed coordinates every 10.0 ps
                                                 ; nstxout-compressed replaces nstxtcout
compressed-x-grps       = System                 ; replaces xtc-grps
; Bond parameters
continuation	        = yes		         ; Restarting after NPT
constraint_algorithm    = lincs	                 ; holonomic constraints
constraints	        = all-bonds	         ; all bonds (even heavy atom-H bonds) constrained
lincs_iter	        = 1		         ; accuracy of LINCS
lincs_order	        = 4		         ; also related to accuracy
; Neighborsearching
cutoff-scheme           = Verlet                 ; REMOVE THIS LINE IF USING GPU
ns_type		        = grid		         ; search neighboring grid cells
nstlist		        = 10	                 ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb	        = 1.0		         ; short-range electrostatic cutoff (in nm)
rvdw		        = 1.0		         ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype	        = PME		         ; Particle Mesh Ewald for long-range electrostatics
pme_order	        = 4		    	 ; cubic interpolation
fourierspacing	        = 0.16		         ; grid spacing for FFT
; Temperature coupling is on
tcoupl		        = V-rescale	         ; modified Berendsen thermostat
tc-grps		        = Protein Non-Protein	 ; two coupling groups - more accurate
tau_t		        = 0.1	  0.1	         ; time constant, in ps
ref_t		        = 300 	  300	         ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl		        = Parrinello-Rahman	 ; Pressure coupling on in NPT
pcoupltype	        = isotropic	         ; uniform scaling of box vectors
tau_p		        = 2.0		         ; time constant, in ps
ref_p		        = 1.0		         ; reference pressure, in bar
compressibility         = 4.5e-5	         ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc		        = xyz		         ; 3-D PBC
; Dispersion correction
DispCorr	        = EnerPres	         ; account for cut-off vdW scheme
; Velocity generation
gen_vel		        = no		         ; Velocity generation is off ";
################################################################################
#                                 SCRIPT                                       #
################################################################################
my $RD  = "\n\033[31m";  # red
my $GR  = "\n\033[32m";  # green
my $YE  = "\n\033[33m";  # yellow
my $BL  = "\n\033[34m";  # blue
my $RS  = "\033[m\n\n";  # reset
################################################################################
#                                FUNCTION                                      #
################################################################################

sub writeout_mdp_files
{
	open(my $emifh, ">", "em.mdp");
	open(my $ionfh, ">", "ions.mdp");
	open(my $nvtfh, ">", "nvt.mdp");
	open(my $nptfh, ">", "npt.mdp");
	open(my $mdrfh, ">", "md.mdp");

	print $emifh $EMI_MDP;
	close $emifh;
	print $ionfh $ION_MDP;
	close $ionfh;
	print $nvtfh $NVT_MDP;
	close $nvtfh;
	print $nptfh $NPT_MDP;
	close $nptfh;
	print $mdrfh $MDR_MDP;
	close $mdrfh;
}

sub check_gromacs_file
{
	if (system("which gmx") == 0) {
		print("$GR GROMACS IS INSTALLED, CONTINUING ... $RS");
	} else {
		die("$RD GROMACS IS NOT INSTALLED ON THIS MACHINE $RS");
	}
}


sub halt_and_prompt
{
	my $input = <STDIN>;
}

my @command_list = (

	$CMD_PDB2GMX,   $CMD_EDITCONF,  $CMD_SOLVATE,     $CMD_GROMPION,
	$CMD_GENION,    $CMD_GROMPPOT,  $CMD_PE_MINIMIZE, $CMD_PLOT_PE,
	$CMD_GROMPNVT,  $CMD_NVT_EQUIB, $CMD_PLOT_NVT,    $CMD_GROMPNPT,
	$CMD_NPT_EQUIB, $CMD_PLOT_NPT,  $CMD_PLOT_DEN,    $CMD_GROMPMD,
	$CMD_SIMULATE,  $CMD_TRJCONV,   $CMD_PLOT_RMSD,   $CMD_PLOT_GYRA,
	$CMD_PLOT_RMSF

);


sub main
{
	check_gromacs_file();

	print("$BL Starting simulation ... $RS");

	print("$BL Writing MDP files ...");

	writeout_mdp_files();

	print(" done $RS");

	print("$GR Starting command execution, hit Ctrl + C to cancel now, or RETURN to continue. $RS");

	halt_and_prompt();

	for my $execute (@command_list)
	{
		system($execute);
		halt_and_prompt();
	}
}

main();

