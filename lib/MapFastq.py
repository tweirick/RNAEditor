#!/usr/bin/python
import argparse
import os
import multiprocessing
from Helper import Helper


class MapFastq(object):
    """ Maps a fastQ file to the given genome/
    """
    def __init__(self, rna_edit_obj):
        """ Constructor """
        self.rnaEdit = rna_edit_obj

        # Set fastq files and check if the qualities need to be converted.
        if self.rnaEdit.params.paired is True:

            fastq_1_is_phred33 = Helper.isPhred33Encoding(
                self.rnaEdit.fastq_files[0],
                1000000,
                self.rnaEdit.logFile,
                self.rnaEdit.textField)

            fastq_2_is_phred33 = Helper.isPhred33Encoding(
                self.rnaEdit.fastq_files[1],
                1000000,
                self.rnaEdit.logFile,
                self.rnaEdit.textField)

            if not fastq_1_is_phred33 or not fastq_2_is_phred33:
                # Convert to phred33.
                self.fastqFile1 = Helper.convertPhred64toPhred33(
                    self.rnaEdit.fastq_files[0],
                    self.rnaEdit.params.output + "_1_phred33.fastq",
                    self.rnaEdit.logFile,
                    self.rnaEdit.textField)

                self.fastqFile2 = Helper.convertPhred64toPhred33(
                    self.rnaEdit.fastq_files[1],
                    self.rnaEdit.params.output + "_2_phred33.fastq",
                    self.rnaEdit.logFile,
                    self.rnaEdit.textField)
            else:
                self.fastqFile1 = self.rnaEdit.fastq_files[0]
                self.fastqFile2 = self.rnaEdit.fastq_files[1]

        elif not self.rnaEdit.params.paired:
            # Convert to phred33.

            fastq_is_phred33 = Helper.isPhred33Encoding(
                self.rnaEdit.fastq_files[0],
                1000000,
                self.rnaEdit.logFile,
                self.rnaEdit.textField)

            if not fastq_is_phred33:
                self.fastqFile1 = Helper.convertPhred64toPhred33(
                    self.rnaEdit.fastq_files[0],
                    self.rnaEdit.params.output + "_1_phred33.fastq",
                    self.rnaEdit.logFile,
                    self.rnaEdit.textField)
            else:
                self.fastqFile = self.rnaEdit.fastq_files[0]

    def printAttributes(self):

        Helper.info("*** MAP READS WITH FOLLOWING ATTRIBUTES ***", self.rnaEdit.logFile, self.rnaEdit.textField)
        if self.rnaEdit.params.paired:
            Helper.info("\t FastQ-File_1: " + self.fastqFile1, self.rnaEdit.logFile, self.rnaEdit.textField)
            Helper.info("\t FastQ-File_2: " + self.fastqFile2, self.rnaEdit.logFile, self.rnaEdit.textField)
        else:
            Helper.info("\t FastQ-File: " + self.fastqFile, self.rnaEdit.logFile, self.rnaEdit.textField)

        Helper.info("\t outfilePrefix:" + self.rnaEdit.params.output, self.rnaEdit.logFile, self.rnaEdit.textField)
        Helper.info("\t refGenome:" + self.rnaEdit.params.refGenome, self.rnaEdit.logFile, self.rnaEdit.textField)
        Helper.info("\t dbsnp:" + self.rnaEdit.params.dbsnp, self.rnaEdit.logFile, self.rnaEdit.textField)
        Helper.info("\t sourceDir:" + self.rnaEdit.params.sourceDir, self.rnaEdit.logFile, self.rnaEdit.textField)
        Helper.info("\t threads:" + self.rnaEdit.params.threads, self.rnaEdit.logFile, self.rnaEdit.textField)
        Helper.info("\t maxDiff:" + self.rnaEdit.params.maxDiff, self.rnaEdit.logFile, self.rnaEdit.textField)
        Helper.info("\t seedDiff:" + self.rnaEdit.params.seedDiff, self.rnaEdit.logFile, self.rnaEdit.textField)
        Helper.info("\t paired:" + str(self.rnaEdit.params.paired), self.rnaEdit.logFile, self.rnaEdit.textField)
        Helper.info("\t keepTemp:" + str(self.rnaEdit.params.keepTemp), self.rnaEdit.logFile, self.rnaEdit.textField)
        Helper.info("\t overwrite:" + str(self.rnaEdit.params.overwrite), self.rnaEdit.logFile, self.rnaEdit.textField)
        Helper.info("", self.rnaEdit.logFile, self.rnaEdit.textField)

    def startAnalysis(self):

        recalibrated_bam_file = self.rnaEdit.params.output+".noDup.realigned.recalibrated.bam"

        if os.path.isfile(recalibrated_bam_file):
            Helper.info(
                "* * * [Skipping] Mapping result File already exists * * *",
                self.rnaEdit.logFile, self.rnaEdit.textField)

            self.rnaEdit.logFile.flush()

            return recalibrated_bam_file

        if self.rnaEdit.params.paired:  # For paired end sequencing.

            # Align first Fastq Reads to the Genome
            sai_file_1 = self.rnaEdit.params.output+"_1.sai"
            cmd = [
                self.rnaEdit.params.sourceDir+"bwa", "aln",
                "-t", self.rnaEdit.params.threads,
                "-n", self.rnaEdit.params.maxDiff,
                "-k", self.rnaEdit.params.seedDiff, self.rnaEdit.params.refGenome, self.fastqFile1
            ]
            Helper.proceedCommand("Align first Reads with BWA", cmd, self.fastqFile1, sai_file_1, self.rnaEdit)
            
            # Align second Fastq Reads to the Genome
            sai_file_2 = self.rnaEdit.params.output+"_2.sai"
            cmd = [
                self.rnaEdit.params.sourceDir+"bwa", "aln",
                "-t", self.rnaEdit.params.threads,
                "-n", self.rnaEdit.params.maxDiff,
                "-k", self.rnaEdit.params.seedDiff,
                self.rnaEdit.params.refGenome, self.fastqFile2
            ]
            Helper.proceedCommand("Align second Reads with BWA", cmd, self.fastqFile2, sai_file_2, self.rnaEdit)
        
            # Convert sai to sam
            sam_file = self.rnaEdit.params.output+".sam"
            cmd = [
                self.rnaEdit.params.sourceDir + "bwa", "sampe",
                "-r", "@RG\tID:bwa\tSM:A\tPL:ILLUMINA\tPU:HiSEQ2000",
                self.rnaEdit.params.refGenome,
                sai_file_1, sai_file_2, self.fastqFile1, self.fastqFile2
            ]
            Helper.proceedCommand("convert sai to sam", cmd, sai_file_1, sam_file, self.rnaEdit)

        else:  # For single end sequencing.
            # Align Fastq Reads to the Genome.
            sai_file = self.rnaEdit.params.output+".sai"
            cmd = [
                self.rnaEdit.params.sourceDir+"bwa", "aln",
                "-t", self.rnaEdit.params.threads,
                "-n", self.rnaEdit.params.maxDiff,
                "-k", self.rnaEdit.params.seedDiff,
                self.rnaEdit.params.refGenome, self.fastqFile
            ]
            Helper.proceedCommand("Align Reads with BWA", cmd, self.fastqFile, sai_file, self.rnaEdit)

            # Convert sai to sam.
            sam_file = self.rnaEdit.params.output+".sam"

            cmd = [
                self.rnaEdit.params.sourceDir + "bwa", "samse",
                "-r", "@RG\tID:bwa\tSM:A\tPL:ILLUMINA\tPU:HiSEQ2000",
                self.rnaEdit.params.refGenome, sai_file, self.fastqFile
            ]

            Helper.proceedCommand("convert sai to sam", cmd, sai_file, sam_file, self.rnaEdit)
        
        # Convert sam to bam.
        unsorted_bam_file = self.rnaEdit.params.output+".unsorted.bam"
        bam_file = self.rnaEdit.params.output+".bam"

        cmd = [self.rnaEdit.params.sourceDir + "samtools", "sort", sam_file, "-o", bam_file]
        Helper.proceedCommand("Sort Bam File", cmd, sam_file, bam_file, self.rnaEdit)

        cmd = [self.rnaEdit.params.sourceDir + "samtools", "index", bam_file]
        Helper.proceedCommand("Index Bam File", cmd, sam_file, bam_file+".bai", self.rnaEdit)
        
        # mark PCR duplicates
        marked_file = self.rnaEdit.params.output+".noDup.bam"
        cmd = [
            "java", "-Xmx16G", "-jar", self.rnaEdit.params.sourceDir + "picard-tools/MarkDuplicates.jar",
            "INPUT=" + bam_file, "OUTPUT=" + marked_file,
            "METRICS_FILE="+self.rnaEdit.params.output+".pcr.metrics",
            "VALIDATION_STRINGENCY=LENIENT", "CREATE_INDEX=true"
        ]
        Helper.proceedCommand("Remove PCR duplicates", cmd, bam_file, marked_file, self.rnaEdit)

        # Identify Target Regions for realignment
        interval_file = self.rnaEdit.params.output+".indels.intervals"
        cmd = [
            "java", "-Xmx16G", "-jar", self.rnaEdit.params.sourceDir + "GATK/GenomeAnalysisTK.jar",
            "-nt", self.rnaEdit.params.threads,
            "-T", "RealignerTargetCreator",
            "-R", self.rnaEdit.params.refGenome,
            "-I", marked_file,
            "-o", interval_file,
            "-l", "ERROR"
        ]
        Helper.proceedCommand("Identify Target Regions for realignment", cmd, bam_file, interval_file, self.rnaEdit)
        
        # Proceed with realignment.
        realigned_file = self.rnaEdit.params.output+".noDup.realigned.bam"
        cmd = [
            "java", "-Xmx16G", "-jar", self.rnaEdit.params.sourceDir + "GATK/GenomeAnalysisTK.jar",
            "-T", "IndelRealigner",
            "-R", self.rnaEdit.params.refGenome,
            "-I", marked_file,
            "-l", "ERROR",
            "-targetIntervals", interval_file,
            "-o", realigned_file
        ]
        Helper.proceedCommand("Proceed Realignment", cmd, interval_file, realigned_file, self.rnaEdit)

        # Find Quality Score recalibration spots.
        recal_file = self.rnaEdit.params.output+".recalSpots.grp"
        cmd = [
            "java", "-Xmx16G", "-jar", self.rnaEdit.params.sourceDir + "GATK/GenomeAnalysisTK.jar",
            "-T", "BaseRecalibrator",
            "-l", "ERROR",
            "-R", self.rnaEdit.params.refGenome,
            "-knownSites", self.rnaEdit.params.dbsnp, "-I", realigned_file,
            "-cov", "CycleCovariate",
            "-cov", "ContextCovariate",
            "-o", recal_file]

        print(" ".join(cmd))

        Helper.proceedCommand("Find Quality Score recalibration spots", cmd, realigned_file, recal_file, self.rnaEdit)
        
        # Proceed quality score recalibration.
        cmd = [
            "java", "-Xmx16G", "-jar", self.rnaEdit.params.sourceDir + "GATK/GenomeAnalysisTK.jar",
            "-T", "PrintReads",
            "-l", "ERROR",
            "-R", self.rnaEdit.params.refGenome,
            "-I", realigned_file,
            "-BQSR", recal_file,
            "-o", recalibrated_bam_file
        ]
        Helper.proceedCommand("Proceed Quality Score recalibration", cmd, recal_file, recalibrated_bam_file, self.rnaEdit)
        
        return recalibrated_bam_file

    def cleanUp(self):
        if not self.rnaEdit.params.keepTemp:
            if os.path.isfile(self.rnaEdit.params.output+".sai"):
                os.remove(self.rnaEdit.params.output+".sai")
            if os.path.isfile(self.rnaEdit.params.output+".sam"):
                os.remove(self.rnaEdit.params.output+".sam")
            if os.path.isfile(self.rnaEdit.params.output+".bam"):
                os.remove(self.rnaEdit.params.output+".bam")
            if os.path.isfile(self.rnaEdit.params.output+".unsorted.bam"):
                os.remove(self.rnaEdit.params.output+".unsorted.bam")
            if os.path.isfile(self.rnaEdit.params.output+".bam.bai"):
                os.remove(self.rnaEdit.params.output+".bam.bai")
            if os.path.isfile(self.rnaEdit.params.output+".indels.intervals"):
                os.remove(self.rnaEdit.params.output+".indels.intervals")
            if os.path.isfile(self.rnaEdit.params.output+".noDup.bam"):
                os.remove(self.rnaEdit.params.output+".noDup.bam")
            if os.path.isfile(self.rnaEdit.params.output+".noDup.bam.bai"):
                os.remove(self.rnaEdit.params.output+".noDup.bam.bai")
            if os.path.isfile(self.rnaEdit.params.output+".noDup.realigned.bam"):
                os.remove(self.rnaEdit.params.output+".noDup.realigned.bam")
            if os.path.isfile(self.rnaEdit.params.output+".noDup.realigned.bai"):
                os.remove(self.rnaEdit.params.output+".noDup.realigned.bai")
            if os.path.isfile(self.rnaEdit.params.output+".recalSpots.grp"):
                os.remove(self.rnaEdit.params.output+".recalSpots.grp")
            # os.remove(self.outfilePrefix+".realigned.marked.recalibrated.bam")


def checkDependencies(args):
    """
    Checks the existence of the necessary packages and tools
    :param args:
    """
    Helper.newline(1)
    Helper.info("CHECK DEPENDENCIES")
    
    # Check if all tools are present.
    if not os.path.isfile(args.sourceDir+"bwa"):
        Helper.error("BWA not found in %s" % args.sourceDir)
    if not os.path.isfile(args.sourceDir+"picard-tools/SortSam.jar"):
        Helper.error("SortSam.jar not found in %s" % args.sourceDir+"picard-tools")
    if not os.path.isfile(args.sourceDir+"picard-tools/MarkDuplicates.jar"):
        Helper.error("MarkDuplicates.jar not found in %s" % args.sourceDir+"picard-tools")
    if not os.path.isfile(args.sourceDir+"GATK/GenomeAnalysisTK.jar"):
        Helper.error("GenomeAnalysisTK.jar not found in %s" % args.sourceDir+"GATK/")
    if not os.path.isfile(args.sourceDir+"samtools"):
        Helper.error("samtools not found in %s" % args.sourceDir)
    if not os.system("java -version") == 0:
        Helper.error("Java could not be found, Please install java")
    
    # Check if all files are present.
    if not os.path.isfile(args.RefGenome):
        Helper.error("Could not find Reference Genome in %s: " % args.RefGenome)
    # Files for BWA
    if not os.path.isfile(args.RefGenome+".amb"):
        Helper.error("Could not find %s.amb" % args.RefGenome)
        Helper.error("run: 'bwa index %s' to create it" % args.RefGenome)
    if not os.path.isfile(args.RefGenome+".ann"):
        Helper.error("Could not find %s.ann" % args.RefGenome)
        Helper.error("run: 'bwa index %s' to create it" % args.RefGenome)
    if not os.path.isfile(args.RefGenome+".bwt"):
        Helper.error("Could not find %s.bwt" % args.RefGenome)
        Helper.error("run: 'bwa index %s' to create it" % args.RefGenome)
    if not os.path.isfile(args.RefGenome+".pac"):
        Helper.error("Could not find %s.pac" % args.RefGenome)
        Helper.error("run: 'bwa index %s' to create it" % args.RefGenome)
    if not os.path.isfile(args.RefGenome+".sa"):
        Helper.error("Could not find %s.sa" % args.RefGenome)
        Helper.error("run: 'bwa index %s' to create it" % args.RefGenome)
    
    # Files for GATK
    if not os.path.isfile(args.RefGenome.replace(".fastq", ".dict")):
        Helper.error("Could not find %s" % args.RefGenome.replace(".fastq", ".dict"))
        Helper.error(
            "run: 'java -jar %s/picard-tools/CreateSequenceDictionary.jar R=%s  O= %s.dict' to create it" %
            (args.sourceDir, args.RefGenome, args.RefGenome.replace(".fastq", ".dict")))
    if not os.path.isfile(args.RefGenome+".fai"):
        Helper.error("Could not find %s.fai" % args.RefGenome)
        Helper.error("run: 'samtools faidx %s' to create it" % args.RefGenome)

    # SNP databases
    if not os.path.isfile(args.dbsnp):
        Helper.error("Could not find %s: " % args.dbsnp)


if __name__ == '__main__':
    # Parse command line arguments and set defaults.
    parser = argparse.ArgumentParser(
        description='map FastQ Files to the given genome and realigns the reads for SNP-calling.')

    parser.add_argument(
        '-i', '--input',
        metavar='Fastq-File',
        type='+',
        help='Input fastq files (maximum two for paire-end-sequencing)',
        required=True
    )
    parser.add_argument(
        "-r", "--RefGenome",
        metavar='Fasta-File',
        help="File that contains the reference sequences",
        type=argparse.FileType('r'),
        default='')
    parser.add_argument(
        '-s', '--dbsnp',
        help=' SNP database (dbSNP) in VCF format (downloaded from the GATK homepage)',
        type=argparse.FileType('r'),
        default='/media/databases/human/dbsnp_135.b37.vcf'
    )
    parser.add_argument(
        '-o', '--output',
        metavar='output-prefix',
        type=str,
        help='prefix that is written in front of the output files', default="default"
    )
    parser.add_argument(
        '-d', '--sourceDir',
        help='- Directory to all the tools [default: /bin/]',
        default='bin/',
        type=Helper.readable_dir
    )
    parser.add_argument(
        '-t', '--threads',
        help='number of threads',
        type=int,
        default=multiprocessing.cpu_count()-1
    )
    parser.add_argument(
        '-n', '--maxDiff',
        help=' maximum Number of mismatches in the reads (int) or error rate in percentage (float)[0.04]',
        type=float, default=0.04
    )
    parser.add_argument(
        '--seedDiff', help='maximum Number of mismatches in the seed sequence (int)[2]',
        type=int,
        default=2
    )
    parser.add_argument(
        '-p', '--paired',
        help="Use this paramater if you have paired end reads [false]",
        action='store_true',
        default=False
    )
    parser.add_argument(
        '--keepTemp',
        help='keep the intermediate Files [False]',
        action='store_true',
        default=False
    )
    parser.add_argument(
        '--overwrite',
        help='overwrite existing Files [False]',
        action='store_true',
        default=False
    )
    
    args = parser.parse_args()
    checkDependencies(args)
    
    map_fastq_obj = MapFastq(
        args.input, args.RefGenome.name, args.dbsnp.name, args.output,
        args.sourceDir, args.threads, args.maxDiff, args.seedDiff,
        args.paired, args.keepTemp, args.overwrite
    )

    map_fastq_obj.startAnalysis()
    del map_fastq_obj
