#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <fstream>
#include <iostream>
#include <string>

int main() {
  std::string input_file = "chembl_33.sdf";
  std::string output_file = "smiles.test";
  std::string smarts = "[#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]";
  RDKit::RWMol *patt = RDKit::SmartsToMol(smarts);
  std::ofstream out;
  out.open(output_file);
  std::unique_ptr<RDKit::ROMol> mol;
  RDKit::SDMolSupplier mol_supplier(input_file, true);
  while (!mol_supplier.atEnd()) {
    mol.reset(mol_supplier.next());

    RDKit::MatchVectType res;
    if (RDKit::SubstructMatch(*mol, *patt, res)) {
      std::cout << "Pattern matched molecule" << std::endl;
    }
    for (size_t i = 0; i < res.size(); ++i) {
      std::cout << "(" << res[i].first << "," << res[i].second << ") ";
    }
    std::cout << std::endl;

    auto smiles = RDKit::MolToSmiles(*mol);
    out << smiles << std::endl;
  }
  out.close();
  return 0;
}
