## Synopsis

Protein-Digester contains the source code for CPPDigester. This program can digest FASTA files containing protein sequences. 
The internal code uses a simple interface to add proteases based on a 8 surrounding amino acids (P1-P4 and P1A-P4A). See the 'proteinases.h' source code. 
Adding your own proteinases should therefore be a trivial task.
The code runs very fast and is fairly simple, but just gets the job done.

## Installation

Requirements:
- CMAKE build automation

Use the available macros in the Makefile or use the cmake itself:

$ mkdir build; cd build; cmake -DCMAKE_BUILD_TYPE=Release ..; make install

## API Reference

The CPPDigester executable is located in the 'bin' folder. Use CPPDigester --help to view the available options and descriptions.

In short:

--minpeplength, -L
	Minimal length of the generated peptides
	
--maxpeplength, -M
	Maximum length of the generated peptides
	
--maxmissed, -C
	Maximum of missed cleavages
	
--proteinase, -P
	Proteinase used for cleaving. Available proteinases: trypsin, arg-c, asp-n, lys-c, chymotrypsin, chymotrypsintrypsin, V8-E, V8-DE
