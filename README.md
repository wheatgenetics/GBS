GBS
===

Tools for GBS Data Management

The programs in the GBS_Utilities folder are used for GBS Sequence Data Management.

GBS Sequence Data Management consists of the following tasks:

1. Download of GBS sequence data files from sequencing facilities.
2. Verify the integrity of downloaded GBS files.
3. Create or rename GBS sequence files to conform to the standard GBS file naming format.
4. Filtering of GBS sequence file to remove short reads of less than 75bp (Nextseq files only).
5. Update wheatgenetics gbs table flowcell, lane, num_lines and md5sum columns for each file.
6. Checking % valid reads in each GBS file and % reads found in each blank well.
7. Checking DNA quantification values for blank wells relative to other wells.
8. Storage of GBS files on Beocat in /bulk/jpoland/sequence GBS sequence file repository
9. Backup of GBS files to external NAS.

Python programs that support these tasks are listed below.

