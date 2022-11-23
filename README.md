# AutoBak
Automated bacteria characterization Workflow

Goal of this repository is to provided an automated workflow 
for Bacteria characterisation.

## Requirements

### Inputs
* a sample sheet linking sample name to Fastq files
* a metadata sheet containing sample information from LIMS (form to define)
* a config file describing execution parameters 
* a QC table with acceptance criteria (species-wise)

### Outputs
* Assemblies
* QC check results
* cgMLST clusters, species-wise
* low priority: characterisation reports (AMR, Virulence, serotypes)
* Species-wise  + global catalog of samples and metadata
* Export-ready tables for external comm (NRLs, EFSA, Land-DB)

### Other reqs
* must be possible to not include samples in databases (comparison mode)
* must be possible to analyze a dataset independently of existing , e.g. for ring trials or validations(independant mode)

## Workflow 
* sanity check inputs
* parse LIMS metadata
* run AQUAMIS on fastq
* compare QC results to acceptance criteria (AutoQC)
* produce an AutoQC report
* create a sample sheet with accepted samples
* run chewieSnake in compare or join mode (to be defined)
* check allel quality and refine sample sheet excluding bad sample
* merge samples in cgMLST databases

## Modules
* LIMS parser
* AQUAMIS wrapper
* AutoQC
* chewieSnake wrapper
* export module
