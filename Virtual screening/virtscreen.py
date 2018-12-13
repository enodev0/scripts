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
                 resfile="ScreeningResults",
                 sampfile="deltaGsamplingData",
                 repeats=2,
                 config="conf.txt",
                 ligand_directory="Ligands/",
                 output_directory="Results/",
                 minimization=False,
                 exhaustiveness=8):

        self.sampfile = sampfile + "-" + str(
            repeats) + "-repeats-" + "-" + str(dt.date.today()) + ".csv"
        self.resfile = resfile + "-" + str(repeats) + "-repeats-" + "-" + str(
            dt.date.today()) + ".csv"
        self.ldir = ligand_directory
        self.rdir = output_directory
        self.minimize = minimization
        self.vinaconf = config
        self.rigor = exhaustiveness
        self.stats = [
        ]  # list of lists to hold mean and SD of each ligand bindng energy
        self.sampling_stats = [
        ]  # list of lists containing binding energy sampling data per ligand
        self.bestbinders = [
        ]  # holds the names of the top five best binding ligands

        header_tmp = ["Ligands", "mean-dG", "stdev-dG", "#samples"]
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

    def convert_ligands(self):
        print("Converting ligands, assuming structures are already in 3D\n\n")
        for file in os.listdir(self.ldir):
            if file.endswith(".sdf"):
                print(
                    "Converting ligand: " + os.path.splitext(file)[0] + "\n\n")
                if self.minimize:
                    print(
                        "Minimizing ligands with steepest descent and GAFF force field ...\n\n"
                    )
                    sp.call(
                        [
                            "obabel -i sdf " + self.ldir + file + " -O " +
                            self.ldir + os.path.splitext(file)[0] + ".sdf" +
                            " --minimize --ff GAFF --sd"
                        ],
                        shell=True)

                else:
                    print("Automatic structure minimization is turned off\n\n")

                sp.call(
                    [
                        "obabel -i sdf " + self.ldir + file + " -O " +
                        self.ldir + os.path.splitext(file)[0] + ".pdbqt" +
                        " -p -partialcharge gasteiger"
                    ],
                    shell=True)

            else:
                print("Could not find ligand SDFs, exiting ...\n\n")
                raise SystemExit()

    def start_screening(self, mode="SINGLE_PASS"):

        print("\nStarting virtual screening run ...\n\n")

        for file in os.listdir(self.ldir):

            if mode == "REPLICATE" and str(
                    os.path.splitext(file)
                [0]) in self.bestbinders[0] and file.endswith(".pdbqt"):

                basename = self.ldir + os.path.splitext(file)[0]
                pathname = self.rdir + os.path.splitext(file)[0]
                os.makedirs(pathname)
                for rep in range(self.repeats):  # from 0 to self.repeats-1
                    os.makedirs(pathname + "/" + str(rep))
                    sp.call(
                        [
                            "vina --config " + self.vinaconf + " --ligand " +
                            self.ldir + file + " --out " + self.rdir +
                            os.path.splitext(file)[0] + "/" + str(rep) + "/" +
                            "out" + "-pass-" + str(rep) + "-" + file +
                            " --log " + self.rdir + os.path.splitext(file)[0] +
                            "/" + str(rep) + "/" + "log-" + "pass-" +
                            str(rep) + "-" + os.path.splitext(file)[0] +
                            ".txt" + " --exhaustiveness " + str(self.rigor)
                        ],
                        shell=True)

            elif mode == "SINGLE_PASS" and file.endswith(".pdbqt"):

                basename = self.ldir + os.path.splitext(file)[0]
                pathname = self.rdir + os.path.splitext(file)[0]
                os.makedirs(pathname)
                for rep in range(self.repeats):  # from 0 to self.repeats-1
                    os.makedirs(pathname + "/" + str(rep))
                    sp.call(
                        [
                            "vina --config " + self.vinaconf + " --ligand " +
                            self.ldir + file + " --out " + self.rdir +
                            os.path.splitext(file)[0] + "/" + str(rep) + "/" +
                            "out" + "-pass-" + str(rep) + "-" + file +
                            " --log " + self.rdir + os.path.splitext(file)[0] +
                            "/" + str(rep) + "/" + "log-" + "pass-" +
                            str(rep) + "-" + os.path.splitext(file)[0] +
                            ".txt" + " --exhaustiveness " + str(self.rigor)
                        ],
                        shell=True)

            else:
                continue

    def process_results(self, mode="SINGLE_PASS"):
        print("\nProcessing results ...\n\n")

        chunks = []  # list containing mean and SD of binding energy / ligand
        rawvars = []  # holds raw binding energies per ligand temporarily
        energies = [
        ]  # list containing raw binding energy values per ligand for t-test

        for file in os.listdir(
                self.ldir):  # file is a builtin, consider changing the name

            if file.endswith(".pdbqt") and mode == "SINGLE_PASS":

                basename = os.path.splitext(file)[0]
                pathname = self.rdir + os.path.splitext(file)[0]
                chunks.append(basename)
                energies.append(basename)

                for rep in range(self.repeats):  # from 0 to self.repeats-1
                    logname = "log-" + "pass-" + str(
                        rep) + "-" + os.path.splitext(file)[0] + ".txt"
                    target = pathname + "/" + str(rep) + "/" + logname
                    # WARNING the below command will fail if vina does not generate 9 models
                    # because it expects the file to have a certain number of lines, which can
                    # vary if the number of models changes. Consider a more robust implementation.
                    cmdline = "tail -n+26 " + target + " | awk '/-/ { print $2 }' | awk 'FNR==2'"

                    res = sp.Popen(
                        [cmdline], shell=True, stdout=sp.PIPE).communicate()[
                            0].strip().decode('ascii')
                    rawvars.append(float(res))
                    energies.append(float(res))

                if self.repeats < 2:
                    # not enough data to calculate mean and stdev
                    # mean will be the single datapoint, stdev will be zero
                    chunks.append(rawvars[0])
                    chunks.append(0)
                    chunks.append(self.repeats)
                else:
                    # we can calculate mean and stdev only if #repeats >= 2
                    chunks.append(stat.mean(rawvars))
                    chunks.append(stat.stdev(rawvars))
                    chunks.append(self.repeats)

                del rawvars[:]
                self.stats.append(chunks.copy())
                del chunks[:]
                self.sampling_stats.append(energies.copy())
                del energies[:]

            elif mode == "REPLICATE" and str(
                    os.path.splitext(file)
                [0]) in self.bestbinders[0] and file.endswith(".pdbqt"):

                basename = os.path.splitext(file)[0]
                pathname = self.rdir + os.path.splitext(file)[0]
                chunks.append(basename)
                energies.append(basename)

                for rep in range(self.repeats):  # from 0 to self.repeats-1
                    logname = "log-" + "pass-" + str(
                        rep) + "-" + os.path.splitext(file)[0] + ".txt"
                    target = pathname + "/" + str(rep) + "/" + logname
                    cmdline = "tail -n+26 " + target + " | awk '/-/ { print $2 }' | awk 'FNR==2'"

                    res = sp.Popen(
                        [cmdline], shell=True, stdout=sp.PIPE).communicate()[
                            0].strip().decode('ascii')
                    rawvars.append(float(res))
                    energies.append(float(res))

                if self.repeats < 2:
                    # not enough data to calculate mean and stdev
                    # mean will be the single datapoint, stdev will be zero
                    chunks.append(rawvars[0])
                    chunks.append(0)
                    chunks.append(self.repeats)
                else:
                    # we can calculate mean and stdev only if #repeats >= 2
                    chunks.append(stat.mean(rawvars))
                    chunks.append(stat.stdev(rawvars))
                    chunks.append(self.repeats)

                del rawvars[:]
                self.stats.append(chunks.copy())
                del chunks[:]
                self.sampling_stats.append(energies.copy())
                del energies[:]

            else:
                continue

    def print_sampling_data(self):
        print("\nPrinting sampling data ...\n\n")

        sg = open(self.sampfile, "w")
        i = 0
        j = len(self.sampling_stats[0])
        for chunk in self.sampling_stats:
            for element in chunk:
                i += 1
                sg.write(str(element))
                if i != j:
                    sg.write(", ")
                else:
                    continue
            i = 0
            sg.write("\n")

        sg.close()

    def print_results(self):
        print("\nPrinting results ...\n\n")

        sr = open(self.resfile, "w")
        i = 0
        j = len(self.stats[0])
        for chunk in self.stats:
            for element in chunk:
                i += 1
                sr.write(str(element))
                if i != j:
                    sr.write(", ")
                else:
                    continue
            i = 0
            sr.write("\n")

        sr.close()

    def extract_bestbinder(self, num=5):
        print("\nSorting results ... \n\n")

        ligand = []  # holds the names of all the ligands screened
        affinity = [
        ]  # holds the corresponding affinities of all the ligands screened

        for i in range(len(self.sampling_stats)):
            ligand.append(self.sampling_stats[i][0])
            affinity.append(self.sampling_stats[i][1])

        affinity, ligand = zip(*sorted(zip(affinity, ligand)))

        self.bestbinders.append(ligand[:num])  # names of top 5 best binders

    def verify_bestbinders(self,
                           outdir="ReplicationResults/",
                           repeats=5,
                           exhaustiveness=20,
                           resfile="ReplicationStats",
                           sampfile="RawReplicationData"):

        print("\nVerifying best binders ...\n\n")

        self.rdir = outdir
        self.repeats = repeats
        self.rigor = exhaustiveness

        self.sampfile = sampfile + "-" + str(
            repeats) + "-repeats-" + "-" + str(dt.date.today()) + ".csv"
        self.resfile = resfile + "-" + str(repeats) + "-repeats-" + "-" + str(
            dt.date.today()) + ".csv"

        del self.stats[:]

        header_tmp = ["Ligands", "mean-dG", "stdev-dG", "#samples"]
        self.stats.append(header_tmp.copy())
        del header_tmp[:]
        del self.sampling_stats[:]

        self.start_screening(mode="REPLICATE")
        self.process_results(mode="REPLICATE")
        self.print_results()
        self.print_sampling_data()

    def split_and_convert_results(self):
        print("Splitting result files ... \n\n")

        for dirs in os.listdir("ReplicationResults/"):
            for dirs_1 in os.listdir("ReplicationResults/" + dirs):
                for files in os.listdir(
                        "ReplicationResults/" + dirs + "/" + dirs_1 + "/"):
                    if files.endswith(".pdbqt"):
                        print("\n\nSplitting file: " + files + "\n\n")
                        target_path = "ReplicationResults/" + dirs + "/" + dirs_1 + "/"
                        sp.call(
                            ["vina_split --input " + target_path + files],
                            shell=True)
                        sp.call(
                            ["rm " + target_path + files],
                            shell=True)  # delete the concatenated file
                        # call openbabel for conversion
                        for files in os.listdir("ReplicationResults/" + dirs +
                                                "/" + dirs_1 + "/"):
                            if files.endswith(".pdbqt"):
                                # extract binding info from PDBQT
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
#    def count_lines(self):
#        file = len([f for f in os.walk("Ligands/").next()[2] if f[-4:] == ".sdf"])
#        print(file)


def main():

    s = Screen(
        repeats=1,
        resfile="screenres-isoform1",
        sampfile="gsampfile-isoform1",
        config="conf.txt",
        ligand_directory="Ligands/",
        output_directory="Results/",
        minimization=False,
        exhaustiveness=10)

    s.convert_ligands()
    s.start_screening()
    s.process_results()
    s.extract_bestbinder(num=5)
    s.print_sampling_data()
    s.print_results()
    s.verify_bestbinders(repeats=5, exhaustiveness=25)
    s.split_and_convert_results()

    raise SystemExit()


if __name__ == "__main__":
    main()
