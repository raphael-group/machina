#include "utils.h"
#include "sampletree.h"
#include "clonetree.h"
#include <fstream>

int main(int argc, char** argv)
{
  if (argc != 8)
  {
    std::cerr << "Usage: " << argv[0] << " <SEED> <#CLONES> <#NRSAMPLES> <S> <lS> <T> <lT> <R>" << std::endl;
    return 1;
  }
  
  int seed = boost::lexical_cast<int>(argv[1]);
  int nrClones = boost::lexical_cast<int>(argv[2]);
  int nrSamples = boost::lexical_cast<int>(argv[3]);
  std::string filenameS = argv[4];
  std::string filenameVertexLabeling = argv[5];
  std::string filenameT = argv[6];
  std::string filenameLeafLabeling = argv[7];
  
  g_rng.seed(seed);

  BoolVector dykeWordS = BinaryTree::randomDyckWord(nrSamples - 1);
  BoolVector dykeWordT = BinaryTree::randomDyckWord(nrClones - 1);
  
  SampleTree S(dykeWordS);
  CloneTree T(dykeWordT, nrSamples);
  
  std::ofstream outS(filenameS.c_str());
  if (!outS.good())
  {
    std::cerr << "Could not open '" << filenameS << "' for reading" << std::endl;
    return 1;
  }
  
  std::ofstream outT(filenameT.c_str());
  if (!outT.good())
  {
    std::cerr << "Could not open '" << filenameT << "' for reading" << std::endl;
    return 1;
  }
  
  std::ofstream outLeafLabeling(filenameLeafLabeling.c_str());
  if (!outLeafLabeling.good())
  {
    std::cerr << "Could not open '" << filenameLeafLabeling << "' for reading" << std::endl;
    return 1;
  }
  
  std::ofstream outVertexLabeling(filenameVertexLabeling.c_str());
  if (!outVertexLabeling.good())
  {
    std::cerr << "Could not open '" << filenameVertexLabeling << "' for reading" << std::endl;
    return 1;
  }
  
  S.writeVertexLabeling(outVertexLabeling);
  outVertexLabeling.close();
  
  S.write(outS);
  outS.close();
  
  T.writeLeafLabeling(outLeafLabeling);
  outLeafLabeling.close();
  
  T.write(outT);
  outT.close();
  
  return 0;
}

int main2(int argc, char** argv)
{
  if (argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " <S> <lS>" << std::endl;
    return 1;
  }
  
  std::string filenameS = argv[1];
  std::ifstream inS(filenameS.c_str());
  if (!inS.good())
  {
    std::cerr << "Could not open '" << filenameS << "' for reading" << std::endl;
    return 1;
  }
  
  std::string filenameVertexLabeling = argv[2];
  std::ifstream inVertexLabeling(filenameVertexLabeling.c_str());
  if (!inVertexLabeling.good())
  {
    std::cerr << "Could not open '" << filenameVertexLabeling << "' for reading" << std::endl;
    return 1;
  }
  
  SampleTree S;
  if (!S.read(inS)) return 1;
  if (!S.readVertexLabeling(inVertexLabeling)) return 1;
  
  std::cout << SampleTree::bracketNotation(S.dyckWord()) << std::endl;
  
  return 0;
}
