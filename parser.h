#ifndef parser_H
#define parser_H

#include <vector>
#include <string>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;

class parser
{
public:
  parser(int ac, char** av) {
    
    status = 0;
    fverbose = false;
    fMaxEvents = -1;
    fMT = false;

    cout << "== Job configuration: " << endl;
    po::options_description desc("Allowed options");

    try {

      //po::options_description desc("Allowed options");
      desc.add_options()
        ("help,h", "produce help message")
        ("job,j", po::value<string>(), "set sample name.")
        ("output-path,o", po::value<string>(), "set path to output directory.")
        ("input-file,i", po::value< vector<string> >()->multitoken(), "set input file(s).")
        ("max-events,m", po::value<int>(), "set a maximum number of events to be run.")
        ("parallel,p", po::value<string>()->implicit_value(""), "enable multi-threading.")
        ("skim,s", po::value<string>(), "write a skimmed tree in a seperate root file.")
        ("verbose,v", po::value<string>()->implicit_value(""), "enable verbosity.")
        ;

      po::variables_map vm;
      po::store(po::parse_command_line(ac, av, desc), vm);
      po::notify(vm);

      
      if (vm.count("help")) {
        cout << "Usage: <executable> [options]\n";
        status = 1;
      }

      if (vm.count("job")) {
        jobname = vm["job"].as<string>();
        cout << "Job name was set to: "
             << jobname << "\n";
      } else {
        cout << "Job name was not set. Use flag -j. Name set to \"test\" \n" << "\n";
        jobname = "test";
      }

      if (vm.count("skim")) {
        skim = vm["skim"].as<string>();
        cout << "Skim enabled, root filename set to: "
             << skim << "\n";
      } 

      if (vm.count("output-path")) {
        outpath = vm["output-path"].as<string>();
        cout << "Output path was set to: "
             << outpath << "\n";
      } else {
        cout << "Output path was not set. Use flag -o. Path set to ./  " << "\n";
      }

      if (vm.count("input-file")) {
        cout << "Input file(s) was set to: " << "\n";
        //<< vm["input-file"].as< vector<string> >() << "\n";

        inputfiles = vm["input-file"].as< vector<string> >();
        for ( vector<string>::iterator it = inputfiles.begin(); it != inputfiles.end(); it++ )
          {
            cout << "  " << *it << endl;
          }

      } else {
        cout << "Input file(s) was not set. This option is required. Use flag -i." << "\n";
        status = 1;
      }

      if (vm.count("verbose")) {
        cout << "Verbosity enabled."
             << "\n";
        fverbose = true;
      }
      
      if (vm.count("parallel")) {
        cout << "Multi-threading enabled."
             << "\n";
        fMT = true;
      }

      if (vm.count("max-events")) {
        fMaxEvents = vm["max-events"].as< int >();
        cout << "Run over maximum "
             << fMaxEvents << " events.\n";
      }

    }
    catch(exception& e) {
      cerr << "== error: " << e.what() << "\n";
      status = 1;
    }
    catch(...) {
      cerr << "== Exception of unknown type!\n";
    }
    
    if (status)
      cout << desc << "\n";

    cout << "==" << endl;

  }
  ~parser() {}
  string JobName() { return jobname; }
  string OutputPath() { return outpath; }
  string SkimName() { return skim; }
  vector<string> InputFiles() { return inputfiles; }
  bool Verbose() { return fverbose; }
  bool MT() { return fMT; }
  bool Exit() { return status; }
private:
  bool status;
  vector<string> inputfiles;
  string outpath;
  string jobname;
  bool fverbose;
  bool fMT;
  int fMaxEvents;
  string skim;
};

#endif
