#include <GraphMol/Atom.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <fstream>
#include <iostream>
#include <string>
#include <limits>

#define DEBUG

std::vector<std::vector<int>>
find_smarts_matches(RDKit::ROMol *rdmol, std::string smarts) {
  RDKit::RWMol *qmol = RDKit::SmartsToMol(smarts);
  // ordered map so we can avoid sorting later
  // this looks exactly like what the python does
  // http://www.rdkit.org/docs/GettingStartedInC%2B%2B.html#atom-map-indices-in-smarts
  std::map<int, unsigned int> idx_map;
  for (auto atom : qmol->atoms()) {
    auto smirks_index = atom->getAtomMapNum();
    if (smirks_index) {
      idx_map[smirks_index - 1] = atom->getIdx();
    }
  }
  std::vector<int> map_list;
  for (const auto &[key, value] : idx_map) {
	std::cout << value << " ";
    map_list.push_back(value);
  }
  std::cout << std::endl;

  std::vector<RDKit::MatchVectType> res;
  RDKit::SubstructMatchParameters params;
  params.useChirality = true;
  params.maxMatches = UINT_MAX;
  params.uniquify = false;
  std::vector<std::vector<int>> ret;
  if (RDKit::SubstructMatch(*rdmol, *qmol, res, &params)) {
    for (size_t i = 0; i < res.size(); ++i) {
      std::vector<int> tmp;
#ifdef DEBUG
      std::cout << "Match " << i + 1 << " : ";
#endif
      for (size_t j = 0; j < map_list.size(); ++j) {
        tmp.push_back(res[i][map_list[j]].second);
#ifdef DEBUG
        std::cout << res[i][map_list[j]].second << " ";
#endif
      }
      ret.push_back(tmp);
#ifdef DEBUG
      std::cout << std::endl;
#endif
    }
  }
#ifdef DEBUG
  std::cout << std::endl;
#endif
  return ret;
}

int main() {
  std::string input_file = "chembl_33.sdf";
  std::string output_file = "smiles.test";
  std::string smarts = "[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]";
  std::ofstream out;
  out.open(output_file);
  std::unique_ptr<RDKit::ROMol> mol;
  RDKit::SDMolSupplier mol_supplier(input_file, true);
  for (size_t i = 0; i < 50 && !mol_supplier.atEnd(); ++i) {
    mol.reset(mol_supplier.next());
    auto smiles = RDKit::MolToSmiles(*mol);
	std::cout << smiles << std::endl;
    find_smarts_matches(&*mol, smarts);

  }
  out.close();
  return 0;
}