# Mummer alignments of Fusarium genomes for Circos plots		

Date: 07/11/22

Nucmer (version 4.0.0rc1) from Mummer was used to align the high quality *cubense*, *lycopersici* and *sacchari* genomes to help identify accessory and core regions, as well as systemic blocks. 

The following Fusarium Assemblies were used in the analysis:

- F._oxysporum_f._sp._cubense_160527.fna
- F._oxysporum_f._sp._cubense_UK0001.fna
- F._oxysporum_f._sp._lycopersici_4287.fna
- F._sacchari_FS66.fna

Comparisons (Each pairwise alignment is stored in the corresponding directory):

- Fol2478_VsFocuUK0001 
- Fol2478_VsFocu160527 

- FocuUK0001_Vs_Foc160527      

- FocuUK0001_Vs_F.SacchariFS66 
- Foc160527_Vs_F.SacchariFS66  


## Nucmer Analysis

The following commands were used for each alignment and filtering. Note this was all run locally using the AlignmentEnv conda environment. 
### Fol2478_VsFocuUK0001

Mummer Commands:

```bash
# Perform analysis
nucmer --maxmatch ../F._oxysporum_f._sp._lycopersici_4287.fna ../F._oxysporum_f._sp._cubense_UK0001.fna 

mv out.delta FoLy_and_FoCubUK0001.delta

delta-filter -g FoLy_and_FoCubUK0001.delta >FoLy_and_FoCubUK0001-deltafilter-g.delta

# plot
mummerplot FoLy_and_FoCubUK0001-deltafilter-g.delta -R ../F._oxysporum_f._sp._lycopersici_4287.fna -Q ../F._oxysporum_f._sp._cubense_UK0001.fna --filter --layout --png

# create circos plot data
show-coords -c -l -L 1000 -r -T FoLy_and_FoCubUK0001-deltafilter-g.delta > FoLy_and_FoCubUK0001-deltafilter-g.coords.txt


# ceate Excel tsv:
cat FoLy_and_FoCubUK0001-deltafilter-g.coords.txt >> Fol4278_UK0001_Links.txt
open -a "Microsoft Excel" Fol4278_UK0001_Links.txt
#The file is then formatted using Excel as a circos links file. 

```

### Fol2478_VsFocu160527

Mummer Commands:

```bash
# Perform analysis
nucmer --maxmatch F._oxysporum_f._sp._lycopersici_4287.fna F._oxysporum_f._sp._cubense_160527.fna

mv out.delta Fol2478_VsFocu160527.delta

delta-filter -g  Fol2478_VsFocu160527.delta > Fol2478_VsFocu160527.deltafilter-g.delta

# plot
mummerplot Fol2478_VsFocu160527.deltafilter-g.delta -R ./F._oxysporum_f._sp._lycopersici_4287.fna -Q F._oxysporum_f._sp._cubense_160527.fna --filter --layout --png

# create circos plot data
show-coords -c -l -L 1000 -r -T Fol2478_VsFocu160527.deltafilter-g.delta > Fol2478_VsFocu160527.deltafilter-g.coords.txt

# Create Excel tsv:
cat Fol2478_VsFocu160527.deltafilter-g.coords.txt >> Fol2478_VsFocu160527_Links.txt
open -a "Microsoft Excel" Fol2478_VsFocu160527_Links.txt
# The file is then formatted using Excel as a circos links file. 
```

### FocuUK0001_Vs_Foc160527

Mummer Commands:

```bash
# perform analysis
nucmer --maxmatch ./F._oxysporum_f._sp._cubense_UK0001.fna ./F._oxysporum_f._sp._cubense_160527.fna

mv out.delta FocuUK0001_Vs_Foc160527.delta

delta-filter -g FocuUK0001_Vs_Foc160527.delta > FocuUK0001_Vs_Foc160527.deltafilter-g.delta

# plot
mummerplot FocuUK0001_Vs_Foc160527.deltafilter-g.delta -R F._oxysporum_f._sp._cubense_UK0001.fna -Q F._oxysporum_f._sp._cubense_160527.fna --filter --layout --png

# create circos plot data
show-coords -c -l -L 1000 -r -T FocuUK0001_Vs_Foc160527.deltafilter-g.delta > FocuUK0001_Vs_Foc160527.deltafilter-g.coords.txt

# Create Excel tsv:
cat FocuUK0001_Vs_Foc160527.deltafilter-g.coords.txt >> FocuUK0001_Vs_Foc160527_Links.txt
open -a "Microsoft Excel" FocuUK0001_Vs_Foc160527_Links.txt
# The file is then formatted using Excel as a circos links file. 
```

### FocuUK0001_Vs_F.SacchariFS66

Mummer Commands:

```bash
# perform analysis
nucmer --maxmatch F._oxysporum_f._sp._cubense_UK0001.fna F._sacchari_FS66.fna 

mv out.delta FocuUK0001_Vs_F.SacchariFS66.delta

delta-filter -g FocuUK0001_Vs_F.SacchariFS66.delta > FocuUK0001_Vs_F.SacchariFS66.deltafilter-g.delta

# plot
mummerplot FocuUK0001_Vs_F.SacchariFS66.deltafilter-g.delta -R F._oxysporum_f._sp._cubense_UK0001.fna -Q F._sacchari_FS66.fna --filter --layout --png

# create circos plot data 
show-coords -c -l -L 1000 -r -T FocuUK0001_Vs_F.SacchariFS66.deltafilter-g.delta > FocuUK0001_Vs_F.SacchariFS66.deltafilter-g.coords.txt

# Create Excel tsv:
cat FocuUK0001_Vs_F.SacchariFS66.deltafilter-g.coords.txt >> FocuUK0001_Vs_F.SacchariFS66_Links.txt
open -a "Microsoft Excel" FocuUK0001_Vs_F.SacchariFS66_Links.txt 
#The file is then formatted using Excel as a circos links file. 
```

### Foc160527_Vs_F.SacchariFS66

Mummer Commands:

```bash
# perfrom analysis
nucmer --maxmatch F._oxysporum_f._sp._cubense_160527.fna F._sacchari_FS66.fna

mv out.delta Foc160527_Vs_F.SacchariFS66.delta

delta-filter -g Foc160527_Vs_F.SacchariFS66.delta > Foc160527_Vs_F.SacchariFS66.deltafilter-g.delta

# plot
mummerplot Foc160527_Vs_F.SacchariFS66.deltafilter-g.delta -R F._oxysporum_f._sp._cubense_160527.fna -Q F._sacchari_FS66.fna --filter --layout --png

# create circos plot data
show-coords -c -l -L 1000 -r -T Foc160527_Vs_F.SacchariFS66.deltafilter-g.delta > Foc160527_Vs_F.SacchariFS66.coords.txt

# Create Excel tsv:
cat Foc160527_Vs_F.SacchariFS66.coords.txt >> Foc160527_Vs_F.SacchariFS66_Links.txt
open -a "Microsoft Excel" 
```

## Circos Analysis

The delta filtered files were then used to build a corresponding links.txt file for circos.

The circos karytopye file was built using the samtools faidx command to create a FASTA index, then formatting into the appropriate circos layout. The Circos karyotype files were then merged from each of the assemblies into a master karytoptye file.