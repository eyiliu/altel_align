// from: Igor Rubinskiy, DESY <mailto:igorrubinsky@gmail.com>
#ifndef EUTELMULTILINEFIT_H
#define EUTELMULTILINEFIT_H

#include <string>
#include <vector>
#include <map>

using IntVec = std::vector<int>;
using FloatVec = std::vector<float>;
using DoubleVec = std::vector<double>;
using StringVec = std::vector<string>;

class Mille;
class EUTelMille{
public:
  EUTelMille ();

  void processEvent (LCEvent * evt);
  void end();


  static void FitTrack(
    unsigned int nPlanesFitter,
    double xPosFitter[],
    double yPosFitter[],
    double zPosFitter[],
    double xResFit[],
    double yResFit[],
    double chi2Fit[2],
    double residXFit[],
    double residYFit[],
    double angleFit[2]
    );


  IntVec _FixedPlanes; //only for internal usage
  IntVec _FixParameter;

  std::string _binaryFilename;
  int _alignMode;

  float _telescopeResolution;

  FloatVec _residualsXMin;
  FloatVec _residualsYMin;
  FloatVec _residualsXMax;
  FloatVec _residualsYMax;

  FloatVec _resolutionX;
  FloatVec _resolutionY;
  FloatVec _resolutionZ;

  std::string _pedeSteerfileName;
  bool _runPede;
  int _usePedeUserStartValues;

  FloatVec _pedeUserStartValuesX;
  FloatVec _pedeUserStartValuesY;
  FloatVec _pedeUserStartValuesZ;
  FloatVec _pedeUserStartValuesAlpha;
  FloatVec _pedeUserStartValuesBeta;
  FloatVec _pedeUserStartValuesGamma;

private:

  //! Event number
  int _iEvt;

  // Statistics
  int _nMilleDataPoints;
  int _nMilleTracks;

  // Mille
  Mille * _mille;

  //! Conversion ID map.
  /*! In the data file, each cluster is tagged with a detector ID
   *  identify the sensor it belongs to. In the geometry
   *  description, there are along with the sensors also "passive"
   *  layers and other stuff. Those are identify by a layerindex. So
   *  we need a conversion table to go from the detectorID to the
   *  layerindex.
   */

  size_t _nPlanes;

  std::vector<DoubleVec > _xPos;
  std::vector<DoubleVec > _yPos;
  std::vector<DoubleVec > _zPos;

  std::vector<DoubleVec > _trackResidX;
  std::vector<DoubleVec > _trackResidY;
  std::vector<DoubleVec > _trackResidZ;

  double * _xPosHere;
  double * _yPosHere;
  double * _zPosHere;
  double * _waferResidX;
  double * _waferResidY;
  double * _waferResidZ;
  double * _telescopeResolX;
  double * _telescopeResolY;
  double * _telescopeResolZ;
  double * _xFitPos;
  double * _yFitPos;

  DoubleVec _siPlaneZPosition;

};
#endif

