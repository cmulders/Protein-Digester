/*
 * fastaProtein.h
 *
 *  Created on: Mar 9, 2015
 *      Author: geforce
 */

#ifndef FASTAPROTEIN_H_
#define FASTAPROTEIN_H_

#include "proteinases.h"

#include <string>

class FastaProtein
{
	private:
		std::string header;
		std::string sequence;

	public:
		FastaProtein(std::string fastaHeader, std::string fastaSequence) :
				header(fastaHeader), sequence(fastaSequence)
		{
		}

		~FastaProtein()
		{
		}

		void digestProtein(std::ostream &outStream, Proteinase *proteinase, unsigned minLength, unsigned maxLength, unsigned maxMissed);


		std::string getHeader()
		{
			return header;
		}

		std::string getSequence()
		{
			return sequence;
		}
};

#endif /* FASTAPROTEIN_H_ */
