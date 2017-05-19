#!/usr/bin/python
"""
Main Class to detect RNA-editing in a given FastQ file
"""
from __future__ import print_function
import sys
from lib.Helper import Helper, Parameters
from lib.createDiagrams import createDiagramms
# from lib import Helper.Helper, Helper.Parameters
#from lib import Helper.Helper as Helper
from lib.MapFastq import MapFastq
from lib.CallEditingSites import CallEditingSites
import multiprocessing, argparse, os
import traceback
import textwrap
import gc
import subprocess


class RnaEdit:
    def __init__(self, fastq_files, params, text_field, out_dir=None, out_base=None):

        if isinstance(params, Parameters):
            self.params = params
        else:
            Helper.error("Params has to be Instance of Parameters")
        try:
            if isinstance(text_field, QtGui.QTextEdit) or text_field == 0:
                self.textField = text_field
            else:
                Helper.error("text_field has to be Instance of QtGui.QTextEdit or 0")
        except NameError:
            self.textField = 0

        self.fastq_files = fastq_files

        # Hold the running Popen object.
        self.runningCommand = False
        self.isTerminated = False
        # Check if the input files are there.

        # Hold basic statistic values of the run.
        basicStatDict = {}

        # Set directory where the outputFiles should be written to
        if self.params.output == "default" and out_dir is None and out_base is None:

            if self.fastq_files[0].endswith("noDup.realigned.recalibrated.bam"):
                self.sampleName = fastq_files[0][
                                  fastq_files[0].rfind("/") +
                                  1:fastq_files[0].rfind(".noDup.realigned.recalibrated.bam")]

                self.out_dir = fastq_files[0][0:fastq_files[0].rfind("/") + 1]
            else:
                self.sampleName = fastq_files[0][fastq_files[0].rfind("/") + 1:fastq_files[0].rfind(".")]

                # out_dir = /path/to/output/rnaEditor/samplename/
                self.out_dir = fastq_files[0][0:fastq_files[0].rfind("/") + 1] + "rnaEditor/" + self.sampleName + "/"
                # An edit to improve command line usability.

        else:
            # Retaining original behavior just in case I missed something.
            if out_base:
                self.sampleName = out_base
            else:
                self.sampleName = fastq_files[0][fastq_files[0].rfind("/") + 1:fastq_files[0].rfind(".")]

            # Retaining original behavior just in case I missed something.
            if out_dir:
                self.out_dir = out_dir
            else:
                self.out_dir = self.params.output.rstrip("/") + "/"

        # Set final parameter path.
        self.params.output = self.out_dir + self.sampleName

        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir, mode=0755)
            os.chmod(self.out_dir, 0755)

        # create folder for html output
        if not os.path.exists(self.out_dir + "/html"):
            os.makedirs(self.out_dir + "/html", mode=0755)
            os.chmod(self.out_dir, 0755)

        self.checkDependencies()

        # Check if the input files are there.
        self.printParameters()

    def run(self):

        #try:
        self.startAnalysis()
        #except Exception:
        #    Helper.error("RnaEditor Failed", self.logFile, self.textField)

        # At this point the RnaEditor has successfully finished.
        cmd = ["python", os.getcwd() + "/createDiagrams.py", "-o", self.params.output]

        a = subprocess.call(cmd)

        self.emit(QtCore.SIGNAL("taskDone"), self.params.output + ".html")

    def startAnalysis(self):
        """
        START MAPPING
        """
        if self.fastq_files[0].endswith("bam"):
            if self.fastq_files[0].endswith("noDup.realigned.recalibrated.bam"):
                Helper.info("Bam File given. Skip mapping", self.logFile, self.textField)
                self.mapFastQ = None
                map_result_file = self.fastq_files[0]
            else:
                Helper.error(
                    "Bam File was not mapped with RnaEditor, this is not supported. " +
                    "Please provide the fastq Files to RnaEditor",
                    self.logFile, self.textField, "red")
        else:
            self.mapFastQ = MapFastq(self)
            map_result_file = self.mapFastQ.startAnalysis()

        """
        START CALLING EDITING SITES
        """
        self.callEditSites = CallEditingSites(map_result_file, self)
        result = self.callEditSites.startAnalysis()

        # Finished.
        self.isTerminated = True

        Helper.status("rnaEditor Finished with %s" % self.params.output,
                      self.logFile, self.textField, "green", True)
        Helper.status("Open %s to see the results" % self.params.output + ".html",
                      self.logFile, self.textField, "green", True)

        self.cleanUp()

        cmd = ["python", os.path.dirname(os.path.realpath(__file__)) + "/createDiagrams.py", "-o", self.params.output]
        a = subprocess.call(cmd)

    def stopSafely(self):
        self.quit()
        Helper.info("Analysis was stopped by User", self.logFile, self.textField)

    def stopImmediately(self):
        if hasattr(self, 'callEditSites'):
            self.callEditSites.cleanUp()
        self.isTerminated = True

        if self.runningCommand != False:
            self.runningCommand.kill()
        else:
            self.terminate()
            self.wait()
        Helper.error("Analysis was terminated by User", self.logFile, self.textField)

    def cleanUp(self):
        # print "deleteAssay " + str(self)
        if self.runningCommand != False:
            self.runningCommand.kill()

        try:
            if self.mapFastQ != None:
                self.mapFastQ.cleanUp()
            del self.mapFastQ
        except AttributeError:
            Helper.error("could not delete MapFastQ instance", self.logFile, self.textField)
        try:
            self.callEditSites.cleanUp()
            del self.callEditSites
        except AttributeError:
            Helper.error("could not delete RnaEdit instance", self.logFile, self.textField)

    def checkDependencies(self):
        """checks if all files are there
        if all programs are installed properly and if the output directory is writable"""
        try:
            self.logFile = open(self.params.output + ".log", "w+")
        except IOError:
            Helper.error("Cannot open Log File", textField=self.textField)

        if type(self.fastq_files) == list:
            self.fastq_files = self.fastq_files
        elif type(self.fastqFile) == str:
            self.fastq_files = [self.fastq_files]
        else:
            Helper.error("FastQ File has wrong variable type", self.logFile, self.textField)

        for file in self.fastq_files:
            if not os.path.isfile(file):
                Helper.error("Could not find: %s" % file, self.logFile, self.textField)

        '''
        Checks the existence of the necessary packages and tools
        :param sourceDir: folder which contains all the software
        '''
        Helper.newline(1)
        Helper.info("CHECK DEPENDENCIES", self.logFile, self.textField)

        # check if all tools are there
        if not os.path.isfile(self.params.sourceDir + "bwa"):
            Helper.error("BWA not found in %s" % self.params.sourceDir, self.logFile, self.textField)
        if not os.path.isfile(self.params.sourceDir + "picard-tools/SortSam.jar"):
            Helper.error("SortSam.jar not found in %s" % self.params.sourceDir + "picard-tools", self.logFile,
                         self.textField)
        if not os.path.isfile(self.params.sourceDir + "picard-tools/MarkDuplicates.jar"):
            Helper.error("MarkDuplicates.jar not found in %s" % self.params.sourceDir + "picard-tools", self.logFile,
                         self.textField)
        if not os.path.isfile(self.params.sourceDir + "GATK/GenomeAnalysisTK.jar"):
            Helper.error("GenomeAnalysisTK.jar not found in %s" % self.params.sourceDir + "GATK/", self.logFile,
                         self.textField)
        if not os.path.isfile(self.params.sourceDir + "blat"):
            Helper.error("blat not found in %s" % self.params.sourceDir, self.logFile, self.textField)
        if not os.path.isfile(self.params.sourceDir + "samtools"):
            Helper.error("samtools not found in %s" % self.params.sourceDir, self.logFile, self.textField)
        if not os.system("java -version") == 0:
            Helper.error("Java could not be found, Please install java", self.logFile, self.textField)

        # check if all files are there
        if not os.path.isfile(self.params.refGenome):
            Helper.error("Could not find Reference Genome in %s: " % self.params.refGenome, self.logFile,
                         self.textField)

        # Files for BWA
        if not os.path.isfile(self.params.refGenome + ".amb"):
            Helper.warning("Could not find %s.amb" % self.params.refGenome, self.logFile, self.textField)
            Helper.error("run: 'bwa index %s' to create it" % self.params.refGenome, self.logFile, self.textField)
        if not os.path.isfile(self.params.refGenome + ".ann"):
            Helper.warning("Could not find %s.ann" % self.params.refGenome, self.logFile, self.textField)
            Helper.error("run: 'bwa index %s' to create it" % self.params.refGenome, self.logFile, self.textField)
        if not os.path.isfile(self.params.refGenome + ".bwt"):
            Helper.warning("Could not find %s.bwt" % self.params.refGenome, self.logFile, self.textField)
            Helper.error("run: 'bwa index %s' to create it" % self.params.refGenome, self.logFile, self.textField)
        if not os.path.isfile(self.params.refGenome + ".pac"):
            Helper.warning("Could not find %s.pac" % self.params.refGenome, self.logFile, self.textField)
            Helper.error("run: 'bwa index %s' to create it" % self.params.refGenome, self.logFile, self.textField)
        if not os.path.isfile(self.params.refGenome + ".sa"):
            Helper.warning("Could not find %s.sa" % self.params.refGenome, self.logFile, self.textField)
            Helper.error("run: 'bwa index %s' to create it" % self.params.refGenome, self.logFile, self.textField)

        # Files for GATK
        if self.params.refGenome.endswith("fasta"):
            if not os.path.isfile(self.params.refGenome.replace(".fasta", ".dict")):
                Helper.warning("Could not find %s" % self.params.refGenome.replace(".fasta", ".dict"), self.logFile,
                               self.textField)
                Helper.error("run: 'java -jar %spicard-tools/CreateSequenceDictionary.jar R=%s  O= %s' to create it" % (
                self.params.sourceDir, self.params.refGenome, self.params.refGenome.replace(".fastq", ".dict")),
                             self.logFile, self.textField)
        elif self.params.refGenome.endswith("fa"):
            if not os.path.isfile(self.params.refGenome.replace(".fa", ".dict")):
                Helper.warning("Could not find %s" % self.params.refGenome.replace(".fa", ".dict"), self.logFile,
                               self.textField)
                Helper.error("run: 'java -jar %spicard-tools/CreateSequenceDictionary.jar R=%s  O= %s' to create it" % (
                self.params.sourceDir, self.params.refGenome, self.params.refGenome.replace(".fa", ".dict")),
                             self.logFile, self.textField)
        else:
            Helper.error("RefGenome has wrong suffix. Either '.fa' or '.fasta'")
        if not os.path.isfile(self.params.refGenome + ".fai"):
            Helper.warning("Could not find %s.sai" % self.params.refGenome, self.logFile, self.textField)
            Helper.error("run: 'samtools faidx %s' to create it" % self.params.refGenome, self.logFile, self.textField)

        # SNP databases
        if not os.path.isfile(self.params.dbsnp):
            Helper.error("Could not find dbSNP database %s: " % self.params.dbsnp, self.logFile, self.textField)
        if not os.path.isfile(self.params.hapmap) and self.params.hapmap != "None":
            Helper.error("Could not find Hapmap database %s: " % self.params.hapmap, self.logFile, self.textField)
        if not os.path.isfile(self.params.omni) and self.params.omni != "None":
            Helper.error("Could not find Omni database %s: " % self.params.omni, self.logFile, self.textField)
        if not os.path.isfile(self.params.esp) and self.params.esp != "None":
            Helper.error("Could not find 1000G database %s: " % self.params.esp, self.logFile, self.textField)

        # region Files
        if not os.path.isfile(self.params.aluRegions):
            Helper.error("Could not find %s: " % self.params.aluRegions, self.logFile, self.textField)

        if not os.path.isfile(self.params.gtfFile):
            Helper.error("Could not find %s: " % self.params.gtfFile, self.logFile, self.textField)

        Helper.info("Dependencies satisfied", self.logFile, self.textField)

        # check if display can activate

    def printParameters(self):

        Helper.info("*** Start RnaEditor with: ***", self.logFile, self.textField)
        if self.fastq_files[0].endswith(".bam"):
            Helper.info("\t Bam File: " + self.fastq_files[0], self.logFile, self.textField)
        else:
            if self.params.paired:
                Helper.info("\t FastQ-File_1: " + self.fastq_files[0], self.logFile, self.textField)
                Helper.info("\t FastQ-File_2: " + self.fastq_files[1], self.logFile, self.textField)
            else:
                Helper.info("\t FastQ-File: " + self.fastq_files[0], self.logFile, self.textField)
        Helper.info("\t outfilePrefix:" + self.params.output, self.logFile, self.textField)
        Helper.info("\t refGenome:" + self.params.refGenome, self.logFile, self.textField)
        Helper.info("\t dbsnp:" + self.params.dbsnp, self.logFile, self.textField)
        Helper.info("\t sourceDir:" + self.params.sourceDir, self.logFile, self.textField)
        Helper.info("\t threads:" + self.params.threads, self.logFile, self.textField)
        Helper.info("\t maxDiff:" + self.params.maxDiff, self.logFile, self.textField)
        Helper.info("\t seedDiff:" + self.params.seedDiff, self.logFile, self.textField)
        Helper.info("\t paired:" + str(self.params.paired), self.logFile, self.textField)
        Helper.info("\t keepTemp:" + str(self.params.keepTemp), self.logFile, self.textField)
        Helper.info("\t overwrite:" + str(self.params.overwrite), self.logFile, self.textField)
        Helper.info("", self.logFile, self.textField)


def main(argv):
    app = QtGui.QApplication(argv)
    mainWindow = GuiView()

    app.setApplicationName("RNAEditor")
    app.setApplicationVersion("0.1")

    app_icon = QtGui.QIcon()
    app_icon.addFile('ui/icons/rnaEditor_16x16.png', QtCore.QSize(16, 16))
    app_icon.addFile('ui/icons/rnaEditor_24x24.png', QtCore.QSize(24, 24))
    app_icon.addFile('ui/icons/rnaEditor_32x32.png', QtCore.QSize(32, 32))
    app_icon.addFile('ui/icons/rnaEditor_48x48.png', QtCore.QSize(48, 48))
    app_icon.addFile('ui/icons/rnaEditor_256x256.png', QtCore.QSize(256, 256))
    app_icon.addFile('ui/icons/rnaEditor_512x512.png', QtCore.QSize(512, 512))
    app_icon.addFile('ui/icons/rnaEditor_1024x1024.png', QtCore.QSize(1024, 1024))

    app.setWindowIcon(app_icon)

    mainWindow.show()

    sys.exit(app.exec_())


if __name__ == '__main__':
    if len(sys.argv) > 1:

        parser = argparse.ArgumentParser(
            prog='RnaEditor',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent('''\
             RnaEditor: easily detect editing sites from deep sequencing data'
            ----------------------------------------------------------------
                run without arguments to start the user interface.
            '''))

        parser.add_argument(
            '-i', '--input',
            metavar='Fastq-Files',
            nargs='+',
            type=str,
            help='Input fastq files (maximum two for paired-end sequencing)',
            required=True
        )

        parser.add_argument(
            '-c', '--conf',
            metavar='Configuration File',
            type=str,
            help='Configuration File used to read Parameters for RnaEditor',
            required=True,
            default='configuration.txt')

        parser.add_argument(
            '-o', '--out',
            metavar='Output Folder',
            type=str,
            required=False,
            default=None,
            help="Location to write RNAEditor's output folder. If not passed the location in the config file will be used."
        )

        parser.add_argument(
            '-b', '--base',
            metavar='Base Name',
            type=str,
            required=False,
            default=None,
            help="The base name for the various output files. Highly recomended for command line use."
        )

        args = parser.parse_args()

        parameters = Parameters(args.conf)

        # Note: args.out, args.base these override default output behavior if passed.

        edit = RnaEdit(args.input, parameters, 0, out_dir=args.out, out_base=args.base)

        edit.run()
        #edit.start()
        #edit.wait()

        del edit

    else:
        from PyQt4 import QtGui, QtCore
        from ui.GuiView import GuiView
        main(sys.argv)
else:
    pass    
