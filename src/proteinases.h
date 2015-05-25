/*
 * proteinases.h
 *
 *  Created on: Mar 3, 2015
 *      Author: geforce
 */

#ifndef PROTEINASES_H_
#define PROTEINASES_H_

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <map>
#include <algorithm>

/*
 * Macro's for better readability
 * Offsets in the sequence for the binding sites locations
 * P4 P3 P2 P1 P1A P2A P3A P4A 	<-- Macro
 * P4 P3 P2 P1 P1' P2' P3' P4' 	<-- Expasy location
 * -3 -2 -1  0  1   2   3   4 	<-- offset in sequence (with check if that sequence exists)
 */
#define ATOFF(OFF) offset >= OFF && seq.at(offset-OFF) //Regular offset
#define ATROFF(ROFF) roffset >= ROFF && seq.at(offset+ROFF) //Reversed offset
#define P4(AA) (ATOFF(3) == AA)
#define P3(AA) (ATOFF(2) == AA)
#define P2(AA) (ATOFF(1) == AA)
#define P1(AA) (ATOFF(0) == AA)

#define P1A(AA) (ATROFF(1) == AA)
#define P2A(AA) (ATROFF(2) == AA)
#define P3A(AA) (ATROFF(3) == AA)
#define P4A(AA) (ATROFF(4) == AA)

class Proteinase
{
	public:
		virtual ~Proteinase()
		{
		}
		;

		//Virtual function that every proteinase will incorporate that has an map of bindings sites
		// Offsets in the sequence for the binding sites locations
		// P4 P3 P2 P1 P1' P2' P3' P4'
		// -3 -2 -1 0  1   2   3   4
		virtual bool isCleavable(std::string const &sequence, const size_t offset, const size_t roffset) = 0;

		//Name of proteinase
		virtual std::string getName() = 0;
		//Rules that this proteinase follows
		virtual std::string getRules() = 0;

		//Determines all the cleavage sites for this protein
		std::vector<std::string::size_type> digest(std::string sequence)
		{
			std::vector<std::string::size_type> cleaveSites;

			//Reverse offset
			size_t roffset = sequence.length() - 1;
			for (size_t offset = 0; offset < sequence.length(); offset++, roffset--)
			{
				if (isCleavable(sequence, offset, roffset))
				{
					//Add cleave site
					cleaveSites.push_back(offset + 1);
				}
			}

			return cleaveSites;
		}

};

/*
 * Enzyme name	P4	P3	P2	P1	P1'	P2'
 * Trypsin 		-	-	-	KR	!P	-
 */
class Trypsin: public Proteinase
{
	public:
		~Trypsin()
		{
		}
		;

		bool isCleavable(std::string const &seq, const size_t offset, const size_t roffset)
		{
			if ((P1('K') || P1('R'))&& !P1A('P')){
				return true;
			}
			return false;
		}

		std::string getName()
		{
			return "trypsin";
		}

		std::string getRules()
		{
			return "http://www.matrixscience.com/help/enzyme_help.html";
		}
};

/*
 * Enzyme name	P4	P3	P2	P1		P1'	P2'
 * Chymotrypsin -	-	-	FYWL	!P	-
 */
class Chymotrypsin: public Proteinase
{
	public:
		~Chymotrypsin()
		{
		}
		;

		bool isCleavable(std::string const &seq, const size_t offset, const size_t roffset)
		{
			if ((P1('F') || P1('Y') || P1('W') || P1('L'))&& !P1A('P')){
				return true;
			}
			return false;
		}

		std::string getName()
		{
			return "chymotrypsin";
		}

		std::string getRules()
		{
			return "http://www.matrixscience.com/help/enzyme_help.html";
		}
};

/*
 * Enzyme name	P4	P3	P2	P1	P1'	P2'
 * Arg-C 		-	-	-	R	!P	-
 */

class ArgC: public Proteinase
{
	public:
		~ArgC()
		{
		}
		;

		bool isCleavable(std::string const &seq, const size_t offset, const size_t roffset)
		{
			if (P1('R') && !P1A('P'))
			{
				return true;
			}
			return false;
		}

		std::string getName()
		{
			return "arg-c";
		}

		std::string getRules()
		{
			return "http://www.matrixscience.com/help/enzyme_help.html";
		}
};
/*
 * Enzyme name	P4	P3	P2	P1	P1'	P2'
 * Asp-N 		-	-	-	-	BD	-
 */

class AspN: public Proteinase
{
	public:
		~AspN()
		{
		}
		;

		bool isCleavable(std::string const &seq, const size_t offset, const size_t roffset)
		{
			if (P1A('D') || P1A('B'))
			{
				return true;
			}
			return false;
		}

		std::string getName()
		{
			return "asp-n";
		}

		std::string getRules()
		{
			return "http://www.matrixscience.com/help/enzyme_help.html";
		}
};

/*
 * Enzyme name	P4	P3	P2	P1	P1'	P2'
 * LysC			-	-	-	K	!P	-
 */
class LysC: public Proteinase
{
	public:
		~LysC()
		{
		}
		;

		bool isCleavable(std::string const &seq, const size_t offset, const size_t roffset)
		{
			if (P1('K') && !P1A('P'))
			{
				return true;
			}
			return false;
		}

		std::string getName()
		{
			return "lys-c";
		}

		std::string getRules()
		{
			return "http://www.matrixscience.com/help/enzyme_help.html";
		}
};

/*
 * Enzyme name			P4	P3	P2	P1		P1'	P2'
 * ChymotrypsinTrypsin	-	-	-	KRFYWL	!P	-
 */
class ChymotrypsinTrypsin: public Proteinase
{
	public:
		~ChymotrypsinTrypsin()
		{
		}
		;

		bool isCleavable(std::string const &seq, const size_t offset, const size_t roffset)
		{
			if ((P1('K') || P1('R') || P1('F') || P1('Y') || P1('W') || P1('L'))&& !P1A('P')){
				return true;
			}
			return false;
		}

		std::string getName()
		{
			return "chymotrypsintrypsin";
		}

		std::string getRules()
		{
			return "http://www.matrixscience.com/help/enzyme_help.html";
		}
};

/*
 * Enzyme name	P4	P3	P2	P1	P1'	P2'
 * GluC (V8-E)	-	-	-	EZ	!P	-
 */
class V8E: public Proteinase
{
	public:
		~V8E()
		{
		}
		;

		bool isCleavable(std::string const &seq, const size_t offset, const size_t roffset)
		{
			if ((P1('E') || P1('Z')) && !P1A('P'))
			{
				return true;
			}
			return false;
		}

		std::string getName()
		{
			return "V8-E";
		}

		std::string getRules()
		{
			return "http://www.matrixscience.com/help/enzyme_help.html";
		}
};

/*
 * Enzyme name	P4	P3	P2	P1		P1'	P2'
 * GluC (V8-DE)	-	-	-	BDEZ	!P	-
 */
class V8DE: public Proteinase
{
	public:
		~V8DE()
		{
		}
		;

		bool isCleavable(std::string const &seq, const size_t offset, const size_t roffset)
		{
			if ((P1('B') || P1('D') || P1('E') || P1('Z')) && !P1A('P'))
			{
				return true;
			}
			return false;
		}

		std::string getName()
		{
			return "V8-DE";
		}

		std::string getRules()
		{
			return "http://www.matrixscience.com/help/enzyme_help.html";
		}
};


class ProteinaseFactory
{
	public:
		static Proteinase* getProteinase(std::string proteinase)
		{
			std::transform(proteinase.begin(), proteinase.end(), proteinase.begin(), ::tolower);
			if (proteinase == "trypsin")
			{
				return new Trypsin();
			}
			else if (proteinase == "chymotrypsintrypsin")
			{
				return new ChymotrypsinTrypsin();
			}
			else if (proteinase == "chymotrypsin")
			{
				return new Chymotrypsin();
			}
			else if (proteinase == "arg-c")
			{
				return new ArgC();
			}
			else if (proteinase == "asp-n")
			{
				return new AspN();
			}
			else if (proteinase == "lys-c")
			{
				return new LysC();
			}
			else if (proteinase == "v8-de")
			{
				return new V8DE();
			}
			else if (proteinase == "v8-e")
			{
				return new V8E();
			}
			else
			{
				throw std::invalid_argument("Proteinase type '" + proteinase + "' not recognized!");
			}
		}
};

#endif
