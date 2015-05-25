/*
 * fastaProtein.cpp
 *
 *  Created on: Mar 9, 2015
 *      Author: geforce
 */

#include "fastaProtein.h"

#include <vector>

void FastaProtein::digestProtein(std::ostream &outStream, Proteinase *proteinase, unsigned minLength, unsigned maxLength, unsigned maxMissed)
{
	std::vector < std::string::size_type > cleavageSites;
	std::string::size_type startCleave, endCleave, cleaveSize;

	cleavageSites = proteinase->digest(this->sequence);

	//Add end and beginning of the sequence
	cleavageSites.insert(cleavageSites.begin(), 0);
	cleavageSites.push_back(this->sequence.length());

	size_t totalCleaveSites = cleavageSites.size();

	for (size_t i = 0; i < (totalCleaveSites - 1); i++)
	{
		//std::cout << i << "\n";
		startCleave = cleavageSites.at(i);
		for (unsigned missed = 0; missed <= maxMissed && (missed + i + 1) < totalCleaveSites; missed++)
		{
			//+1 because we want the current cleave position to the next (=1) + missed cleave sites
			endCleave = cleavageSites.at(i + missed + 1);
			cleaveSize = endCleave - startCleave;
			if (cleaveSize >= minLength && (maxLength == 0 || cleaveSize <= maxLength))
			{
				std::string peptide(this->sequence.substr(startCleave, cleaveSize));
				/*
				 * Check if this sequence contains any non-amino acids
				 * 21st and 22nd amino acids	3-Letter	1-Letter
				 * Selenocysteine				Sec			U
				 * Pyrrolysine					Pyl			O
				 *
				 *	In addition to the specific amino acid codes, placeholders are used in
				 *	cases where chemical or crystallographic analysis of a peptide or protein
				 *	cannot conclusively determine the identity of a residue.
				 *	Ambiguous Amino Acids				3-Letter	1-Letter
				 *	Asparagine or aspartic acid			Asx			B
				 *	Glutamine or glutamic acid			Glx			Z
				 *	Leucine or Isoleucine				Xle			J
				 *	Unspecified or unknown amino acid	Xaa			X
				 */

				if (peptide.find_first_of("UOBZJX") == std::string::npos)
				{
					outStream << peptide << "\n";
				}
			}
		}
	}
}

