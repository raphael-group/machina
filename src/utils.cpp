#include "utils.h"
#include <fstream>

std::mt19937 g_rng(0);

std::istream& getline(std::istream& is, std::string& t)
{
  // copied from: http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
  t.clear();
  
  // The characters in the stream are read one-by-one using a std::streambuf.
  // That is faster than reading them one-by-one using the std::istream.
  // Code that uses streambuf this way must be guarded by a sentry object.
  // The sentry object performs various tasks,
  // such as thread synchronization and updating the stream state.
  
  std::istream::sentry se(is, true);
  std::streambuf* sb = is.rdbuf();
  
  for(;;) {
    int c = sb->sbumpc();
    switch (c) {
      case '\n':
        return is;
      case '\r':
        if(sb->sgetc() == '\n')
          sb->sbumpc();
        return is;
      case EOF:
        // Also handle the case when the last line has no line ending
        if(t.empty())
          is.setstate(std::ios::eofbit);
        return is;
      default:
        t += (char)c;
    }
  }
}

bool parseMigrationGraph(const std::string& migrationGraphFile,
                         const StringSet& samples,
                         StringPairList& forcedComigrations)
{
  std::ifstream inFile(migrationGraphFile);
  if (!inFile.good())
  {
    std::cerr << "Could not open '" << migrationGraphFile << "' for reading" << std::endl;
    return false;
  }
  else
  {
    int idx = 0;
    while (inFile.good())
    {
      ++idx;
      std::string line;
      getline(inFile, line);
      
      if (line.empty()) continue;
      
      StringVector s;
      boost::split(s, line, boost::is_any_of("\t "));
      
      if (s.size() < 2)
      {
        std::cerr << "Line " << idx << " is invalid in '" << migrationGraphFile
        << "'" << std::endl;
        return false;
      }
      if (samples.count(s[0]) != 1)
      {
        std::cerr << "Line " << idx << " is invalid in '" << migrationGraphFile
        << "'. Sample '" << s[0] << "' is incorrect." << std::endl;
        return false;
      }
      if (samples.count(s[1]) != 1)
      {
        std::cerr << "Line " << idx << " is invalid in '" << migrationGraphFile
        << "'. Sample '" << s[1] << "' is incorrect." << std::endl;
        return false;
      }
      forcedComigrations.push_back(StringPair(s[0], s[1]));
    }
  }
  
  return true;
}
