/*
 * fastaReader.cpp
 *
 *  Created on: Mar 5, 2015
 *      Author: geforce
 */

#include "fastaReader.h"

//Go the next entry
bool FastaReader::next()
{
	//Delete the contents of the previous entry
	this->currentHeader.clear();
	this->currentSequence.clear();

	if (this->fastaFileStream.eof())
	{
		//Reached end of file, return false
		return false;
	}

	std::string lineBuffer;

	//Read till fasta header ('>') or end of file
	while (!this->fastaFileStream.eof() && std::getline(this->fastaFileStream, lineBuffer) && lineBuffer[0] != '>')
		;

	//End of file
	if (this->fastaFileStream.eof())
	{
		//Reached end of file, return false, no sequence data
		return false;
	}

	this->currentHeader.append(lineBuffer);

	//Read sequence till end of file or next line starts with fasta header ('>')
	while (!this->fastaFileStream.eof() && this->fastaFileStream.peek() != '>' && std::getline(this->fastaFileStream, lineBuffer))
	{
		//Append the fasta sequence to our entry
		this->currentSequence.append(lineBuffer);
	}

	if (this->currentHeader.empty() || this->currentSequence.empty()) return false;

	return true;
}
