#!/usr/bin/env python


"""
Automated virtual screening & analysis with Autodock Vina
"""
"""
REQUIREMENTS:

Python 3.6.4 :: Anaconda, Inc.
AutoDock Vina 1.1.2 (May 11, 2011)
Open Babel 2.3.2 -- Dec 18 2015
awk, tail, cat


Copyright (C) 2018, Somdeb Chatterjee
LICENSE: MIT
"""

import datetime as dt
import statistics as stat
import subprocess as sp
import os


class Screen:
    """
    The screen class represents a virtual screening session
    """

    def __init__(self,
                 resfile="affinitites_with_stats",
                 sampfile="affinities_only",
                 repeats=1,
                 config="conf.txt",
                 ligand_directory="Ligands/",
                 output_directory="Results/",
                 repeats_directory="Replications/",
                 minimization=False,
                 exhaustiveness=8):

        self.sampfile = sampfile + "-" + str(
            repeats) + "-repeats-" + "-" + str(dt.date.today()) + ".csv"
        self.resfile = resfile + "-" + str(repeats) + "-repeats-" + "-" + str(
            dt.date.today()) + ".csv"
        self.ligand_dir = ligand_directory
        self.result_dir = output_directory
        self.repdir = repeats_directory
        self.minimize = minimization
        self.vinaconf = config
        self.rigor = exhaustiveness
        self.stats = [
        ]  # list of lists to hold mean and SD of each ligand bindng energy
        self.sampling_stats = [
        ]  # list of lists containing binding energy sampling data per ligand
        self.bestbinders = [
        ]  # holds the names of the top five best binding ligands
        self.energies = [
        ]  # list containing raw binding energy values per ligand for t-test

        self.chunks = [
        ]  # list containing mean and SD of binding energy / ligand
        self.rawvars = []  # holds raw binding energies per ligand temporarily

        header_tmp = ["Ligand", "mean-dG", "stdev-dG", "#samples"]
        self.stats.append(header_tmp.copy())
        del header_tmp[:]

        if repeats >= 1:
            self.repeats = repeats
        else:
            print("Setting default repeat count (=1)\n\n")
            self.repeats = 1

        if (os.system("which obabel") == 0 and os.system("which vina") == 0):
            print(
                "AutoDock Vina and OpenBabel are available, continuing ...\n\n"
            )
        else:
            print(
                "\n\nFATAL: Check that AutoDock Vina and OpenBabel are installed\n\n"
            )
            raise SystemExit()


    def preprocess_ligands(self):
        print("Converting ligands, assuming structures are already in 3D\n\n")
        for filename in os.listdir(self.ligand_dir):
            if filename.endswith(".sdf"):
                print("Converting ligand: " + os.path.splitext(filename)[0] +
                      "\n\n")
                if self.minimize:
                    print(
                        "Minimizing ligands with steepest descent and GAFF force field ...\n\n"
                    )
                    sp.call(
                        [
                            "obabel -i sdf " + self.ligand_dir + filename +
                            " -O " + self.ligand_dir +
                            os.path.splitext(filename)[0] + ".sdf" +
                            " --minimize --ff GAFF --sd"
                        ],
                        shell=True)

                else:
                    print("Automatic structure minimization is turned off\n\n")

                sp.call(
                    [
                        "obabel -i sdf " + self.ligand_dir + filename + " -O "
                        + self.ligand_dir + os.path.splitext(filename)[0] +
                        ".pdbqt" + " -p -partialcharge gasteiger"
                    ],
                    shell=True)

            else:
                print("Could not find ligand SDFs, exiting ...\n\n")
                raise SystemExit()


    def dock_molecule(self, ligand_name):

        pathname = self.result_dir + os.path.splitext(ligand_name)[0]
        os.makedirs(pathname)

        for rep in range(self.repeats):  # from 0 to self.repeats-1
            os.makedirs(pathname + "/" + str(rep))
            sp.call(
                [
                    "vina --config " + self.vinaconf + " --ligand " +
                    self.ligand_dir + ligand_name + " --out " + pathname +
                    "/" + str(rep) + "/" + "out" + "-pass-" + str(rep) + "-" +
                    ligand_name + " --log " + pathname + "/" + str(rep) + "/" +
                    "log_" + "pass_" + str(rep) + "_" +
                    os.path.splitext(ligand_name)[0] + ".txt" +
                    " --exhaustiveness " + str(self.rigor)
                ],
                shell=True)


    def start_virtual_screen(self, mode="SINGLE_PASS"):

        print("\nStarting virtual screening run ...\n\n")

        for filename in os.listdir(self.ligand_dir):

            if mode == "REPLICATED" and str(
                    os.path.splitext(filename)
                [0]) in self.bestbinders[0] and filename.endswith(".pdbqt"):

                self.dock_molecule(filename)

            elif mode == "SINGLE_PASS" and filename.endswith(".pdbqt"):

                self.dock_molecule(filename)

            else:
                continue


    def extract_data_from_pdbqt(self, filename):

        basename = os.path.splitext(filename)[0]
        pathname = self.result_dir + os.path.splitext(filename)[0]
        self.chunks.append(basename)
        self.energies.append(basename)

        for rep in range(self.repeats):  # from 0 to self.repeats-1
            logname = "log_" + "pass_" + str(rep) + "_" + os.path.splitext(
                filename)[0] + ".txt"
            target = pathname + "/" + str(rep) + "/" + logname
            # WARNING: the below command will fail if vina does not generate 9 models
            # because it expects the file to have a certain number of lines, which can
            # vary if the number of models changes. Consider a more robust implementation.
            # Extract the binding energy from a Vina log file.
            cmdline = "tail -n+26 " + target + " | awk '/-/ { print $2 }' | awk 'FNR==1'"

            res = sp.Popen(
                [cmdline], shell=True,
                stdout=sp.PIPE).communicate()[0].strip().decode('ascii')
            self.rawvars.append(float(res))
            self.energies.append(float(res))

        if self.repeats < 2:
            # not enough data to calculate mean and stdev
            # mean will be the single datapoint, stdev will be zero
            self.chunks.append(self.rawvars[0])
            self.chunks.append(0)
            self.chunks.append(self.repeats)
        else:
            # we can calculate mean and stdev only if #repeats >= 2
            self.chunks.append(stat.mean(self.rawvars))
            self.chunks.append(stat.stdev(self.rawvars))
            self.chunks.append(self.repeats)

        del self.rawvars[:]
        self.stats.append(self.chunks.copy())
        del self.chunks[:]
        self.sampling_stats.append(self.energies.copy())
        del self.energies[:]


    def extract_binding_affinities(self, mode="SINGLE_PASS"):
        print("\nExtracting affinities ...\n\n")

        for filename in os.listdir(self.ligand_dir):

            if filename.endswith(".pdbqt") and mode == "SINGLE_PASS":

                self.extract_data_from_pdbqt(filename)

            elif mode == "REPLICATED" and str(
                    os.path.splitext(filename)
                [0]) in self.bestbinders[0] and filename.endswith(".pdbqt"):

                self.extract_data_from_pdbqt(filename)

            else:
                continue


    def print_results(self, mode="SIMPLE"):
        print("\nPrinting results ...\n\n")

        if mode == "SIMPLE":
            outfile = self.sampfile
            datasrc = self.sampling_stats
        elif mode == "WITH_STATS":
            outfile = self.resfile
            datasrc = self.stats
        else:
            print(
                "\nInvalid mode: "+mode+
                " in print_results()\n\n"
                )
            raise SystemExit()

        filename = open(outfile, "w")
        i = 0
        j = len(datasrc[0])
        for chunk in datasrc:
            for element in chunk:
                i += 1
                filename.write(str(element))
                if i != j:
                    filename.write(", ")
                else:
                    continue
            i = 0
            filename.write("\n")

        filename.close()


    def determine_bestbinders(self, num=5):
        print("\nSorting results ... \n\n")

        ligand = []  # holds the names of all the ligands screened
        affinity = [
        ]  # holds the corresponding affinities of all the ligands screened

        for i in range(len(self.sampling_stats)):
            ligand.append(self.sampling_stats[i][0])
            affinity.append(self.sampling_stats[i][1])

        affinity, ligand = zip(*sorted(zip(affinity, ligand)))

        self.bestbinders.append(ligand[:num])  # names of top 5 best binders


    def verify_bestbinders(
            self,
            repeats=5,
            exhaustiveness=20,
            resfile="replication_affinities_with_stats",
            sampfile="replication_affinities_only"):

        print("\nVerifying best binders ...\n\n")

        self.result_dir = self.repdir
        self.repeats = repeats
        self.rigor = exhaustiveness

        self.sampfile = sampfile + "_" + str(repeats) + "_repeats_" + str(
            dt.date.today()) + ".csv"
        self.resfile = resfile + "_" + str(repeats) + "_repeats_" + str(
            dt.date.today()) + ".csv"

        del self.stats[:]

        header_tmp = ["Ligands", "mean-dG", "stdev-dG", "#samples"]
        self.stats.append(header_tmp.copy())
        del header_tmp[:]
        del self.sampling_stats[:]

        self.start_virtual_screen(mode="REPLICATED")
        self.extract_binding_affinities(mode="REPLICATED")
        self.print_results(mode="SIMPLE")
        self.print_results(mode="WITH_STATS")


    def split_and_convert_pdbqts(self, mode="SINGLE_PASS"):
        print("\nSplitting result files and (PDBQT -->  SDF) ... \n\n")

        target_dir = ""

        if mode == "SINGLE_PASS":
            target_dir = self.result_dir
        elif mode == "REPLICATED":
            target_dir = self.repdir
        else:
            print("\nInvalid mode: "+mode+" in split_and_convert_pdbqts()\n\n")
            raise SystemExit()

        for dirs in os.listdir(target_dir):
            for dirs_1 in os.listdir(target_dir + dirs):
                for files in os.listdir(
                        target_dir + dirs + "/" + dirs_1 + "/"):
                    if files.endswith(".pdbqt"):
                        print("\n\nSplitting file: " + files + "\n\n")
                        target_path = target_dir + dirs + "/" + dirs_1 + "/"
                        sp.call(
                            ["vina_split --input " + target_path + files],
                            shell=True)
                        sp.call(
                            ["rm " + target_path + files],
                            shell=True)  # delete the concatenated file
                        # call openbabel for conversion
                        for files in os.listdir(
                                target_dir + dirs + "/" + dirs_1 + "/"):
                            if files.endswith(".pdbqt"):
                                # extract binding energy, RMS min and max from PDBQT
                                dG = sp.Popen(
                                    [
                                        "cat " + target_path + files +
                                        "| awk '/-/ { print $4 }' | awk 'NR==1{print $1}'"
                                    ],
                                    shell=True,
                                    stdout=sp.PIPE).communicate()[
                                        0].strip().decode('ascii')
                                rms_min = sp.Popen(
                                    [
                                        "cat " + target_path + files +
                                        "| awk '/-/ { print $5 }' | awk 'NR==1{print $1}'"
                                    ],
                                    shell=True,
                                    stdout=sp.PIPE).communicate()[
                                        0].strip().decode('ascii')
                                rms_max = sp.Popen(
                                    [
                                        "cat " + target_path + files +
                                        "| awk '/-/ { print $6 }' | awk 'NR==1{print $1}'"
                                    ],
                                    shell=True,
                                    stdout=sp.PIPE).communicate()[
                                        0].strip().decode('ascii')

                                sp.call(
                                    [
                                        "obabel -i pdbqt " + target_path +
                                        files + " -O " + target_path + dirs +
                                        "-affinity=" + dG + "_rmsMAX=" +
                                        rms_max + "_rmsMIN=" + rms_min +
                                        "_Split.sdf"
                                    ],
                                    shell=True)


def main():

    s = Screen(
        resfile="affinities_with_stats",
        sampfile="affinities",
        config="conf.txt",
        ligand_directory="Ligands/",
        output_directory="Results/",
        repeats_directory="Replications/",
        minimization=False,
        exhaustiveness=8)

    s.preprocess_ligands()
    s.start_virtual_screen()
    s.extract_binding_affinities()
    s.determine_bestbinders(num=3)

    s.print_results(mode="SIMPLE")
    s.print_results(mode="WITH_STATS")
    s.split_and_convert_pdbqts(mode="SINGLE_PASS")

    s.verify_bestbinders(repeats=3, exhaustiveness=8)
    s.split_and_convert_pdbqts(mode="REPLICATED")

    print("\n\n\nAll done!\n\n")
    raise SystemExit()


if __name__ == "__main__":
    main()
