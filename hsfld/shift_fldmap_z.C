#include <TSystem.h>

#include <fstream>
#include <iomanip>
#include <iostream>

namespace {
constexpr const char *kDefaultInput =
    "hsfld/param_history/ShsFieldMap_20210526_Extrapolated";
constexpr const char *kDefaultOutput =
    "hsfld/param_history/ShsFieldMap_20210526_Extrapolated_HS0";
}

void shift_fldmap_z(const char *inputPath = kDefaultInput,
                    const char *outputPath = kDefaultOutput,
                    double zShift = 171.95)
{
  std::ifstream fin(inputPath);
  if (!fin) {
    std::cerr << "Cannot open input file: " << inputPath << std::endl;
    return;
  }

  std::ofstream fout(outputPath);
  if (!fout) {
    std::cerr << "Cannot open output file: " << outputPath << std::endl;
    return;
  }

  int nx = 0;
  int ny = 0;
  int nz = 0;
  double x0 = 0.0;
  double y0 = 0.0;
  double z0 = 0.0;
  double dx = 0.0;
  double dy = 0.0;
  double dz = 0.0;
  fin >> nx >> ny >> nz >> x0 >> y0 >> z0 >> dx >> dy >> dz;

  if (!fin) {
    std::cerr << "Bad header in input file: " << inputPath << std::endl;
    return;
  }

  fout << std::setprecision(12);
  fout << nx << " " << ny << " " << nz << " " << x0 << " " << y0 << " "
       << z0 + zShift << " " << dx << " " << dy << " " << dz << "\n";

  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double bx = 0.0;
  double by = 0.0;
  double bz = 0.0;
  Long64_t nrow = 0;

  while (fin >> x >> y >> z >> bx >> by >> bz) {
    fout << x << " " << y << " " << z + zShift << " " << bx << " " << by
         << " " << bz << "\n";
    ++nrow;
  }

  if (!fin.eof()) {
    std::cerr << "Stopped before EOF while reading: " << inputPath
              << std::endl;
    return;
  }

  std::cout << "Wrote shifted field map: " << outputPath << std::endl;
  std::cout << "Rows shifted: " << nrow << std::endl;
  std::cout << "Header z0: " << z0 << " -> " << z0 + zShift << std::endl;
  std::cout << "Applied z shift: " << zShift << " cm" << std::endl;
}
