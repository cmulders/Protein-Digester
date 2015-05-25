#include "proteinases.h"
#include "fastaReader.h"

#define VERSION_MAJOR 0
#define VERSION_MINOR 2

#define PROGRESSBAR_MAX 40

#include <algorithm>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/exception/diagnostic_information.hpp>

// A helper function to simplify outputing vectors to string streams.
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
	copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
	return os;
}

int main(int argc, char **argv)
{
	namespace po = boost::program_options;

	//Map that holds our command line options
	po::variables_map vm;

	try
	{
		std::string appName = boost::filesystem::basename(argv[0]);

		// Declare a group of options that will be
		// allowed only on command line
		po::options_description generic("Generic options");
		generic.add_options() //
		("version,v", "print version string")("help,h", "produce help message");

		// Declare a group of options that will be
		// allowed both on command line and in
		// config file
		po::options_description config("Configuration");
		config.add_options() //
		("minpeplength,L", po::value<unsigned>()->default_value(1), "Minimal length of the generated peptides") //
		("maxpeplength,M", po::value<unsigned>()->default_value(0), "Maximum length of the generated peptides") //
		("maxmissed,C", po::value<unsigned>()->default_value(0), "Maximum of missed cleavages") //
		("proteinase,P", po::value<std::string>()->default_value("trypsin"),
				"Proteinase used for cleaving.\n" "Available proteinases: trypsin, arg-c, asp-n, lys-c, chymotrypsin, chymotrypsintrypsin, V8-E, V8-DE"); //

		// Hidden options, will be allowed both on command line and
		// in config file, but will not be shown to the user.
		po::options_description hidden("Hidden options");
		hidden.add_options() //
		("input-file", po::value<std::vector<std::string> >(), "input file");

		//Create cmd_line options group for parsing
		po::options_description cmdline_options;
		cmdline_options.add(generic).add(config).add(hidden);

		//Create visible options (without input-file) for display
		po::options_description visible("Allowed options");
		visible.add(generic).add(config);

		po::positional_options_description p;
		p.add("input-file", -1);

		//Load the command line arguments in our variable map and show any errors
		try
		{
			po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
		} catch (po::unknown_option &e)
		{
			std::cout << e.what() << "\n";
			std::cout << visible << "\n";
			return EXIT_FAILURE;
		} catch (std::exception &e)
		{
			std::cout << e.what() << "\n";
			return EXIT_FAILURE;
		}
		//Show help if asked
		if (vm.count("help"))
		{
			std::cout << "Usage: " << appName << " [options] file1.fasta ..." << "\n";
			std::cout << visible << "\n";
			return EXIT_SUCCESS;
		}

		//Print out the version
		if (vm.count("version"))
		{
			std::cout << "CPP digester, version " << VERSION_MAJOR << "." << VERSION_MINOR << "\n";
			return EXIT_SUCCESS;
		}

		//Print out the input files
		if (vm.count("input-file"))
		{
			bool valid = true;

			//Check all the files if they are all valid
			for (auto &file : vm["input-file"].as<std::vector<std::string> >())
			{
				if (boost::filesystem::extension(file) != ".fasta")
				{
					std::cerr << file << " has no vaild extension (" << boost::filesystem::extension(file)
							<< " != .fasta)" << std::endl;
					valid = false;
				}
			}

			//Can print this vector because we overloaded the << for vector printing (see above the main() )
			std::cout << "Input files (n=" << vm["input-file"].as<std::vector<std::string> >().size() << ") are:\n"
					<< vm["input-file"].as<std::vector<std::string> >() << "\n\n";

			if (!valid) return EXIT_FAILURE;
		}
		else
		{
			//If not input files specified
			std::cout << "No input files.\n";
			return EXIT_SUCCESS;
		}

	} catch (std::exception &e)
	{
		std::cout << e.what() << "\n";
		return EXIT_FAILURE;
	}

	/*
	 * Options are loaded, try to recognize our proteinase
	 */

	Proteinase* proteinase;
	try
	{
		//Try to create the proteinase, however we can fail and get a exception back
		proteinase = ProteinaseFactory::getProteinase(vm["proteinase"].as<std::string>());
	} catch (std::invalid_argument &e)
	{
		//Did not have that proteinase in our factory method
		std::cout << e.what() << "\n";
		std::cout << "Use --help for help." << std::endl;
		return EXIT_FAILURE;
	} catch (std::exception &e)
	{
		//Other exception
		std::cout << e.what() << "\n";
		return EXIT_FAILURE;
	}

	std::vector<std::string> inputFiles = vm["input-file"].as<std::vector<std::string> >();

	//Inform the user of the used proteinase
	std::cout << "Starting digestion with '" << proteinase->getName() << std::endl;

	FastaReader *curFastaFile;
	for (auto it = inputFiles.begin(); it != inputFiles.end(); it++)
	{
		std::string outFile(
				boost::filesystem::path(*it).stem().string() + "." + proteinase->getName() + ".peptides.txt");
		//Open the output file and delete current contents
		std::ofstream outStream(outFile, std::ofstream::out | std::ofstream::trunc);

		std::cout << "Processing: " << (*it) << " outputting in: " << outFile << std::endl;

		//All user output should be indented hereafter
		try
		{
			curFastaFile = new FastaReader(*it);
		} catch (std::exception &e)
		{
			curFastaFile = NULL;
			std::cerr << "\tCould not process current fasta file. Because this error:\n";
			std::cerr << "\t" << e.what() << "\n";
			std::cerr << "\tSkipping this file" << std::endl;
			continue;
		}
		std::cout << "\tFasta filesize: " << curFastaFile->fastaFileSize << " bytes " //
				<< "(" << (curFastaFile->fastaFileSize / 1000 / 1000) << "MB)" << std::endl;

		//For progress bar
		int _last;

		std::vector<std::string::size_type> cleaveSites;
		FastaProtein* curProtein;
		while (curFastaFile->next())
		{
			//Digest the current protein with the proteinase and output to the file stream
			curFastaFile->getProtein()->digestProtein(outStream, proteinase, vm["minpeplength"].as<unsigned>(), vm["maxpeplength"].as<unsigned>(),
					vm["maxmissed"].as<unsigned>());

			//Update progress bar
			int cur(std::ceil(curFastaFile->getProgress() * PROGRESSBAR_MAX));
			if (_last != cur)
			{
				_last = cur;
				std::cerr << std::fixed << std::setprecision(2) //
						<< "\r\t[" << std::string(cur, '#') << std::string(PROGRESSBAR_MAX - cur, ' ') << "] " //
						<< 100 * curFastaFile->getProgress() << "%"; //

			}
		}
		//Update progress bar to 100%
		std::cout << std::fixed << std::setprecision(2) //
				<< "\r\t[" << std::string(PROGRESSBAR_MAX, '#') << "] " << 100.0f << "%" << std::endl;

		std::cout << "\tDigested filesize: " << outStream.tellp() << " bytes (" << (outStream.tellp() / 1000 / 1000)
				<< "MB)" << std::endl;

		std::cout << std::endl;
		delete curFastaFile;
		outStream.close();
	}

	std::cout << "\n" << "Done processing. Exiting" << std::endl;
	return EXIT_SUCCESS;
}
