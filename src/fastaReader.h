/*
 * fastaReader.h
 *
 *  Created on: Mar 5, 2015
 *      Author: geforce
 */

#ifndef FASTAREADER_H_
#define FASTAREADER_H_

#include "fastaProtein.h"

#include <fstream>
#include <stdexcept>



class FastaReader
{
	private:
		std::ifstream fastaFileStream;
		std::string currentHeader;
		std::string currentSequence;

		void setFileSize()
		{
			if (this->fastaFileStream.is_open())
			{
				//Get current pointer in the file
				std::ifstream::pos_type currentOffset = this->fastaFileStream.tellg();

				//Go to the end and record the size
				this->fastaFileStream.seekg(0, std::ios::end);
				this->fastaFileSize = this->fastaFileStream.tellg();

				//Return to the original position
				this->fastaFileStream.seekg(currentOffset, std::ios::beg);
			}
		}


		public:
		std::size_t fastaFileSize;

		FastaReader(std::string fastaFileName) :
				fastaFileSize(0), currentHeader(""), currentSequence("")
		{
			this->fastaFileStream.open(fastaFileName);
			if (!this->fastaFileStream.is_open())
			{
				throw std::runtime_error("Could not open file '" + fastaFileName + "'");
			}
			this->setFileSize();

		}

		~FastaReader()
		{
			if (this->fastaFileStream.is_open())
			{
				this->fastaFileStream.close();
			}
		}

		bool next();

		float getProgress()
		{
			if(this->fastaFileStream.eof())
			{
				return 1.0f;
			}
			return (float) this->fastaFileStream.tellg() / this->fastaFileSize;
		}

		FastaProtein* getProtein()
		{
			return new FastaProtein(this->currentHeader,this->currentSequence);
		}

};

#endif /* FASTAREADER_H_ */
