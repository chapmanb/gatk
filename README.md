This is a small fork of the Broad's Genome Analysis Toolkit (GATK) framework
to provide a ready to run shell script executable as `gatk-framework`. This is
the MIT licensed, publicly available code, so does not include Variant Calling,
Base Quality Score Recalibration, Variant Quality Score Recalibration or other
GATK functionality that requires a license from Broad or Appistry.

The shell script, `gatk-framework`, uses exactly the same parameters as the
standard GATK jar and also correctly passes java specific arguments like `-Xmx`
`-Xms` and `-D` to java.

Releases, coinciding with GATK releases, are available from:
https://github.com/chapmanb/gatk/releases
They contain the shell wrapper and pre-built uberjar with
GATK and all dependencies.

Original code and functionality from Broad: http://www.broadinstitute.org/gatk/
