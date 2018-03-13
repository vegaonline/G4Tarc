#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iterator>
#include <unordered_map>

std::string inFileName = "neutSpec.dat";
std::string outFileName = "neutronET.dat";
std::ifstream ndataIN(inFileName, std::ios::in);
std::ofstream ndataOUT(outFileName, std::ios::out);
std::vector<double> energyTab;
std::vector<double> timeTab;
std::vector<int> countTab;
std::vector <double> :: iterator it;
std::vector <double> :: iterator ie;
std::vector <int> :: iterator ic;
double errLim = 0.1;



template<typename T1, typename T2, typename Allocator>
std::unordered_map< T1, std::size_t > frequencies( std::vector<T1, Allocator> const& src, T2 key, T2 errlimit ) {
  int count = 0;
  std::unordered_map< T1, std::size_t > retval;
  for (auto&& x:src){
    if ( std::abs(x - key) <= errlimit)
      ++retval[x];
  }
  return retval;
}


int main() {
  std::string dummy;
  double myTime, myEnergy;


  std::cout << " Reading the data from the " << inFileName <<  std::endl;
  while (ndataIN){
    getline(ndataIN, dummy);
    std::stringstream ssdummy(dummy);
    ssdummy >> myTime >> myEnergy;
    energyTab.push_back(myEnergy);
    timeTab.push_back(myTime);
  }
  ndataIN.close();

  // choose key energy
  double keyEnergy = 1.0;  // eV
  std::cout << " Total data read <energy values>: " << energyTab.size() << " <time values>: " << timeTab.size() << std::endl;
  std::unordered_map<double, std::size_t> thisMap = frequencies(energyTab, keyEnergy, errLim);
  for (const auto& n: thisMap){
    countTab.push_back(n.second);
  }

  for (it = timeTab.begin(), ie = energyTab.begin(), ic=countTab.begin();
       it != timeTab.end();
       it++, ie++, ic++){
         if ( std::abs(*ie - keyEnergy) <= errLim) {
           std::cout << *it << "  " << *ie << "  " << *ic << std::endl;
         }
    //std::cout << *it << "  " << *ic << std::endl;
  }
  ndataOUT.close();

  return 0;
}
