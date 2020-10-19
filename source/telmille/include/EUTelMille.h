// from: Igor Rubinskiy, DESY <mailto:igorrubinsky@gmail.com>
#ifndef EUTELMULTILINEFIT_H
#define EUTELMULTILINEFIT_H

#include <string>
#include <vector>
#include <map>
#include "myrapidjson.h"

using IntVec = std::vector<int>;
using FloatVec = std::vector<float>;
using DoubleVec = std::vector<double>;
using StringVec = std::vector<std::string>;

class Mille;
class EUTelMille{
public:
  EUTelMille();
  ~EUTelMille();

  void setGeometry(const JsonValue& js);
  void setResolution(double resolX, double resolY);

  void startMilleBinary(const std::string& path);
  void endMilleBinary();

  void fillTrackXYRz(const JsonValue& js);
  void createPedeStreeringModeXYRz(const std::string& path);

  // static void FitTrack(
  //   unsigned int nPlanesFitter,
  //   double xPosFitter[],
  //   double yPosFitter[],
  //   double zPosFitter[],
  //   double xResFit[],
  //   double yResFit[],
  //   double chi2Fit[2],
  //   double residXFit[],
  //   double residYFit[],
  //   double angleFit[2]
  //   );

  static void FitTrack(unsigned int nMeasures,
                       const std::vector<double>& xPosMeasure,
                       const std::vector<double>& yPosMeasure,
                       const std::vector<double>& zPosMeasure,
                       const std::vector<double>& xResolMeasure,
                       const std::vector<double>& yResolMeasure,
                       double& xOriginLine,
                       double& yOriginLine,
                       double& xAngleLine,
                       double& yAngleLine,
                       double& xChisqLine,
                       double& yChisqLine,
                       std::vector<double>& residXFit,
                       std::vector<double>& residYFit);


  // static void FitTrack(unsigned int nPlanesFit,
  //                      const std::vector<double>& xPosFitter,
  //                      const std::vector<double>& yPosFitter,
  //                      const std::vector<double>& zPosFitter,
  //                      const std::vector<double>& xResolFitter,
  //                      const std::vector<double>& yResolFitter,
  //                      std::vector<double>& chi2Fit,
  //                      std::vector<double>& residXFit,
  //                      std::vector<double>& residYFit,
  //                      std::vector<double>& angleFit);

private:
  std::unique_ptr<Mille> m_mille;
  std::string m_binPath;

  size_t m_nPlanes;
  std::map<size_t, size_t> m_indexDet;
  double m_xResolution;
  double m_yResolution;

  std::map<size_t, double> m_xPosDet;
  std::map<size_t, double> m_yPosDet;
  std::map<size_t, double> m_zPosDet;

  std::map<size_t, double> m_alphaPosDet;
  std::map<size_t, double> m_betaPosDet;
  std::map<size_t, double> m_gammaPosDet;
};
#endif

