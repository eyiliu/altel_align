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


namespace gbl{
   class GblDetectorLayer;
}

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


  static std::unique_ptr<gbl::GblDetectorLayer> CreateLayerSit(const std::string& aName, unsigned int layer,
                                                               double xPos, double yPos,
                                                               double zPos, double thickness,
                                                               double uAngle, double uRes,
                                                               double vAngle, double vRes);

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

  std::map<size_t, std::unique_ptr<gbl::GblDetectorLayer>> m_dets;

};




#endif

