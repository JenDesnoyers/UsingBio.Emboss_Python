# UsingBio.Emboss_Python
This script takes a GenBank query search, extracts six nucleotide sequence results and then uses the sequences in a series of Bio.Emboss applications such as: FNeighborCommandline, FDNADistCommandline, NeedleallCommandline, and IepCommandline. 

////////////////////////////////////////////////////////////////////////Notes//////////////////////////////////////////////////////////////////////////////////

**Important: if Emboss is not installed to "usr/local/bin" then you will have to find where each application is on your computer and change these paths where each application is called in the script. Also, it would be fairly easy to add new applications to this script, but be aware that not all applications work. 

When you run the script, it will ask you for what query you would like to use so there is no need to add that to the script itself. 

//////////////////////////////////////////////////////////////////////Dependencies////////////////////////////////////////////////////////////////////////////// 
#Python 3 
#BioPython 
#Bio.Emoboss -> Emboss can be downloaded from here: http://emboss.sourceforge.net/download/ 
