#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <fstream>
#include <iostream>
#include <string>

int main() {
  std::string input_file = "chembl_33.sdf";
  std::string output_file = "smiles.out";
  std::ofstream out;
  out.open(output_file);
  std::unique_ptr<RDKit::ROMol> mol;
  RDKit::SDMolSupplier mol_supplier(input_file, true);
  while (!mol_supplier.atEnd()) {
    mol.reset(mol_supplier.next());
    auto smiles = RDKit::MolToSmiles(*mol);
    out << smiles << std::endl;
  }
  out.close();
  return 0;
}
