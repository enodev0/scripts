#!/usr/bin/env perl

# Automatically dump PDBs at fixed intervals of time from a GROMACS simulation trajectory using gmx trjconv

# 22 -> SELECT GROUP: non-water, includes ions+ligands+proteins

use strict;
use warnings;

my $dump_timestep_pico = 10000;    #picoseconds

my $current_time = 0;

my $timestep_nano = 0;

for ( my $i = 0 ; $i <= 10 ; $i++ ) {
        $timestep_nano = $current_time / 1000;    # ps2ns
        system(
"{ echo 22; echo 22; } | gmx trjconv -s md.tpr -f mdnopbc.xtc -n index.ndx -o complex-at-$timestep_nano-ns.pdb -center -ur compact -pbc mol -dump $current_time"
        );

        print("\n\nDumped $timestep_nano ns PDB.\n\n");
        sleep(3); # seconds
        $current_time = $current_time + $dump_timestep_pico;
}
