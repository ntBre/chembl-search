#include <GraphMol/Atom.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

using namespace std;
using namespace RDKit;
using namespace RDKit::MolOps;

vector<vector<int>> find_smarts_matches(RWMol *rdmol, string smarts) {
  RWMol *qmol = SmartsToMol(smarts);
  cout << MolToSmiles(*qmol) << endl;
  // ordered map so we can avoid sorting later
  // this looks exactly like what the python does
  // http://www.rdkit.org/docs/GettingStartedInC%2B%2B.html#atom-map-indices-in-smarts
  map<int, unsigned int> idx_map;
  for (auto atom : qmol->atoms()) {
    auto smirks_index = atom->getAtomMapNum();
    if (smirks_index) {
      idx_map[smirks_index - 1] = atom->getIdx();
    }
  }
  vector<int> map_list;
  cout << "map list:" << endl;
  for (const auto &[key, value] : idx_map) {
    cout << value << " ";
    map_list.push_back(value);
  }
  cout << endl;

  vector<MatchVectType> res;
  SubstructMatchParameters params;
  params.useChirality = true;
  params.maxMatches = UINT_MAX;
  params.uniquify = false;
  vector<vector<int>> ret;
  if (SubstructMatch(*rdmol, *qmol, res)) {
    for (size_t i = 0; i < res.size(); ++i) {
      vector<int> tmp;
      cout << "Match " << i + 1 << " : ";
      for (size_t j = 0; j < map_list.size(); ++j) {
        tmp.push_back(res[i][map_list[j]].second);
        cout << res[i][map_list[j]].second << " ";
      }
      ret.push_back(tmp);
      cout << endl;
    }
  } else {
    cout << "No matches" << endl;
  }
  cout << endl;
  return ret;
}

int main() {
  string smiles = "Cc1cc(-c2csc(N=C(N)N)n2)cn1C";
  string smarts = "[*:1]-[#16X2,#16X3+1:2]-[#6:3]~[*:4]";
  string input_file = "chembl_33.sdf";
  SDMolSupplier mol_supplier(input_file, true);
  // ROMol *mol = SmilesToMol(smiles);
  ROMol *mol = mol_supplier.next();
  RDKit::RWMol *mol4(new RDKit::RWMol(*mol));
  unsigned int failed;
  sanitizeMol(*mol4, failed,
              SANITIZE_ALL ^ SANITIZE_ADJUSTHS ^ SANITIZE_SETAROMATICITY);
  cout << failed << endl;
  setAromaticity(*mol4, AROMATICITY_MDL);
  assignStereochemistry(*mol4);
  addHs(*mol4);
  cout << MolToSmiles(*mol) << endl;
  find_smarts_matches(mol4, smarts);
  return 0;
}
