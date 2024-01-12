#include <GraphMol/Atom.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <fstream>
#include <iostream>
#include <string>

// #define DEBUG

std::vector<std::tuple<int, int, int, int>>
find_smarts_matches(RDKit::ROMol *rdmol, std::string smarts) {
  RDKit::RWMol *qmol = RDKit::SmartsToMol(smarts);
  // ordered map so we can avoid sorting later
  // this looks exactly like what the python does
  // http://www.rdkit.org/docs/GettingStartedInC%2B%2B.html#atom-map-indices-in-smarts
  std::map<int, int> idx_map;
  for (auto atom : qmol->atoms()) {
    auto smirks_index = atom->getAtomMapNum();
    if (smirks_index != 0) {
      idx_map[smirks_index - 1] = atom->getIdx();
    }
  }
  std::vector<int> map_list;
  for (const auto &[key, value] : idx_map) {
    map_list.push_back(value);
  }

  std::vector<RDKit::MatchVectType> res;
  bool useChirality = true;
  std::vector<std::tuple<int, int, int, int>> ret;
  if (RDKit::SubstructMatch(*rdmol, *qmol, res, useChirality)) {
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
      ret.push_back(std::make_tuple(tmp[0], tmp[1], tmp[2], tmp[3]));
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
    find_smarts_matches(&*mol, smarts);

    auto smiles = RDKit::MolToSmiles(*mol);
    out << smiles << std::endl;
  }
  out.close();
  return 0;
}
