// from: Igor Rubinskiy, DESY <mailto:igorrubinsky@gmail.com>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <numeric>
#include <iostream>
#include <fstream>

#include <cstring>
#include <cmath>


#include "streamlog/streamlog.h"
#include "pstream.h"
#include "Mille.h"
#include "EUTelMille.h"

using namespace std;



/*! Performs analytic straight line fit.
 *
 * Determines parameters of a straight line passing through the 
 * measurement points.
 *
 * @see http://www.desy.de/~blobel/eBuch.pdf page 162
 *
 */
void EUTelMille::FitTrack(unsigned int nPlanesFitter,
                          double xPosFitter[], double yPosFitter[], double zPosFitter[],
                          double xResFitter[], double yResFitter[], double chi2Fit[2],
                          double residXFit[], double residYFit[], double angleFit[2]) {

  int sizearray;
  sizearray = nPlanesFitter;

  double * xPosFit = new double[sizearray];
  double * yPosFit = new double[sizearray];
  double * zPosFit = new double[sizearray];
  double * xResFit = new double[sizearray];
  double * yResFit = new double[sizearray];

  int nPlanesFit = 0;

  for (unsigned int n = 0; n < nPlanesFitter; n++) {


    xPosFit[nPlanesFit] = xPosFitter[n];
    yPosFit[nPlanesFit] = yPosFitter[n];
    zPosFit[nPlanesFit] = zPosFitter[n];
    xResFit[nPlanesFit] = xResFitter[n];
    yResFit[nPlanesFit] = yResFitter[n];
    nPlanesFit++;
  }

  int counter;

  float S1[2]   = {0,0};
  float Sx[2]   = {0,0};
  float Xbar[2] = {0,0};

  float * Zbar_X = new float[nPlanesFit];
  float * Zbar_Y = new float[nPlanesFit];
  for (counter = 0; counter < nPlanesFit; counter++){
    Zbar_X[counter] = 0.;
    Zbar_Y[counter] = 0.;
  }

  float Sy[2]     = {0,0};
  float Ybar[2]   = {0,0};
  float Sxybar[2] = {0,0};
  float Sxxbar[2] = {0,0};
  float A2[2]     = {0,0};

  // define S1
  for( counter = 0; counter < nPlanesFit; counter++ ){
    S1[0] = S1[0] + 1/pow(xResFit[counter],2);
    S1[1] = S1[1] + 1/pow(yResFit[counter],2);
  }

  // define Sx
  for( counter = 0; counter < nPlanesFit; counter++ ){
    Sx[0] = Sx[0] + zPosFit[counter]/pow(xResFit[counter],2);
    Sx[1] = Sx[1] + zPosFit[counter]/pow(yResFit[counter],2);
  }

  // define Xbar
  Xbar[0]=Sx[0]/S1[0];
  Xbar[1]=Sx[1]/S1[1];

  // coordinate transformation !! -> bar
  for( counter = 0; counter < nPlanesFit; counter++ ){
    Zbar_X[counter] = zPosFit[counter]-Xbar[0];
    Zbar_Y[counter] = zPosFit[counter]-Xbar[1];
  }

  // define Sy
  for( counter = 0; counter < nPlanesFit; counter++ ){
    Sy[0] = Sy[0] + xPosFit[counter]/pow(xResFit[counter],2);
    Sy[1] = Sy[1] + yPosFit[counter]/pow(yResFit[counter],2);
  }

  // define Ybar
  Ybar[0]=Sy[0]/S1[0];
  Ybar[1]=Sy[1]/S1[1];

  // define Sxybar
  for( counter = 0; counter < nPlanesFit; counter++ ){
    Sxybar[0] = Sxybar[0] + Zbar_X[counter] * xPosFit[counter]/pow(xResFit[counter],2);
    Sxybar[1] = Sxybar[1] + Zbar_Y[counter] * yPosFit[counter]/pow(yResFit[counter],2);
  }

  // define Sxxbar
  for( counter = 0; counter < nPlanesFit; counter++ ){
    Sxxbar[0] = Sxxbar[0] + Zbar_X[counter] * Zbar_X[counter]/pow(xResFit[counter],2);
    Sxxbar[1] = Sxxbar[1] + Zbar_Y[counter] * Zbar_Y[counter]/pow(yResFit[counter],2);
  }

  // define A2
  A2[0]=Sxybar[0]/Sxxbar[0];
  A2[1]=Sxybar[1]/Sxxbar[1];

  // Calculate chi sqaured
  // Chi^2 for X and Y coordinate for hits in all planes
  for( counter = 0; counter < nPlanesFit; counter++ ){
    chi2Fit[0] += pow(-zPosFit[counter]*A2[0]
                      +xPosFit[counter]-Ybar[0]+Xbar[0]*A2[0],2)/pow(xResFit[counter],2);
    chi2Fit[1] += pow(-zPosFit[counter]*A2[1]
                      +yPosFit[counter]-Ybar[1]+Xbar[1]*A2[1],2)/pow(yResFit[counter],2);
  }

  for( counter = 0; counter < static_cast< int >(nPlanesFitter); counter++ ) {
    residXFit[counter] = (Ybar[0]-Xbar[0]*A2[0]+zPosFitter[counter]*A2[0])-xPosFitter[counter];
    residYFit[counter] = (Ybar[1]-Xbar[1]*A2[1]+zPosFitter[counter]*A2[1])-yPosFitter[counter];
  }

  // define angle
  angleFit[0] = atan(A2[0]);
  angleFit[1] = atan(A2[1]);

  // clean up
  delete [] zPosFit;
  delete [] yPosFit;
  delete [] xPosFit;
  delete [] yResFit;
  delete [] xResFit;

  delete [] Zbar_X;
  delete [] Zbar_Y;

}

void EUTelMille::setGeometry(const JsonValue& js) {
    if(!js.HasMember("geometry")){
    std::fprintf(stderr, "unable to find \"geomerty\" key for detector geomerty from JS\n");
    return;
  }

  const auto &js_geo = js["geometry"];
  const auto &js_dets = js_geo["detectors"];
  std::map<double, size_t> zmap_sort;
  for(const auto& js_det: js_dets.GetArray()){
    size_t id = js_det["id"].GetUint();
    double cz = js_det["center"]["z"].GetDouble();
    zmap_sort[cz]=id;
  }

  m_nPlanes=zmap_sort.size();
  size_t n=0;
  for(auto [the_cz, the_id]: zmap_sort){
    for(const auto& js_det: js_dets.GetArray()){
      size_t id = js_det["id"].GetUint();
      if(id == the_id){
        double cx = js_det["center"]["x"].GetDouble();
        double cy = js_det["center"]["y"].GetDouble();
        double cz = js_det["center"]["z"].GetDouble();
        double rx = js_det["rotation"]["x"].GetDouble();
        double ry = js_det["rotation"]["y"].GetDouble();
        double rz = js_det["rotation"]["z"].GetDouble();
        m_xPosDet[id]=cx;
        m_yPosDet[id]=cy;
        m_zPosDet[id]=cz;

        m_alphaPosDet[id]=rx;
        m_betaPosDet[id]=ry;
        m_gammaPosDet[id]=rz;
        m_indexDet[id] = n;
        n++;
      }
    }
  }
}


void EUTelMille::startMilleBinary(const std::string& path){
  m_mille.reset(new Mille(path.c_str()));
  m_binPath=path;
}


void EUTelMille::endMilleBinary(){
  m_mille.reset();
}

void EUTelMille::fillTrackXYRz(const JsonValue& js) {
  std::vector<double> xPosHit(m_nPlanes,0);
  std::vector<double> yPosHit(m_nPlanes,0);
  std::vector<double> zPosHit(m_nPlanes,0);

  //TODO: set
  std::vector<double> xResolHit(m_nPlanes, 0.02);
  std::vector<double> yResolHit(m_nPlanes, 0.02);

  std::vector<size_t> idHit(m_nPlanes,0);

  if(js.Size()!= m_nPlanes){
    std::fprintf(stderr, "hits number is less than detector number \n");
    throw;
  }
  for(const auto& js_hit : js.GetArray()){
    size_t id =js_hit["id"].GetUint();
    size_t detN = m_indexDet.at(id);
    double xPosDet = m_xPosDet.at(id);
    double yPosDet = m_yPosDet.at(id);
    double zPosDet = m_zPosDet.at(id);

    xPosHit[detN]=xPosDet+js_hit["x"].GetDouble();
    yPosHit[detN]=yPosDet+js_hit["y"].GetDouble();
    zPosHit[detN]=zPosDet;
    idHit[detN]=id;
  }

  std::vector<double> chisqTrack(2,0);
  std::vector<double> angleTrack(2,0);
  std::vector<double> xResidHit(m_nPlanes,0);
  std::vector<double> yResidHit(m_nPlanes,0);

  // Calculate residuals
  FitTrack(m_nPlanes,
           xPosHit.data(),
           yPosHit.data(),
           zPosHit.data(),
           xResolHit.data(),
           yResolHit.data(),
           chisqTrack.data(),
           xResidHit.data(),
           yResidHit.data(),
           angleTrack.data());

  const int nLC = 4; // number of local parameters
  const int nGL = m_nPlanes * 3; // number of global parameters

  std::vector<float> derLC(nLC, 0);
  std::vector<float> derGL(nGL, 0);
  std::vector<int> label(nGL, 0);
  std::iota(std::begin(label), std::end(label), 1);

  // loop over all planes
  for (unsigned int n = 0; n < m_nPlanes; n++) {
    float residual;
    float sigma;

    derGL[((n * 3) + 0)] = -1;
    derGL[((n * 3) + 2)] = yPosHit[n];
    derLC[0] = 1;
    derLC[2] = zPosHit[n];
    residual = xResidHit[n];
    sigma    = xResolHit[n];
    m_mille->mille(nLC,derLC.data(),nGL,derGL.data(),label.data(),residual,sigma);

    derGL[((n * 3) + 0)] = 0;
    derGL[((n * 3) + 2)] = 0;
    derLC[0] = 0;
    derLC[2] = 0;

    derGL[((n * 3) + 1)] = -1;
    derGL[((n * 3) + 2)] = -1 * xPosHit[n];
    derLC[1] = 1;
    derLC[3] = zPosHit[n];
    residual = yResidHit[n];
    sigma    = yResolHit[n];
    m_mille->mille(nLC,derLC.data(),nGL,derGL.data(),label.data(),residual,sigma);

    derGL[((n * 3) + 1)] = 0;
    derGL[((n * 3) + 2)] = 0;
    derLC[1] = 0;
    derLC[3] = 0;

  } // end loop over all planes

  m_mille->end();
}


void EUTelMille::createPedeStreeringModeXYRz(const std::string& path){

  streamlog_out ( MESSAGE4 ) << endl << "Generating the steering file for the pede program..." << endl;

  ofstream steerFile;
  steerFile.open(path.c_str());

  steerFile << "Cfiles" << endl;
  steerFile << m_binPath << endl;
  steerFile << endl;

  steerFile << "Parameter" << endl;

  int counter = 0;
  // loop over all planes
  for (unsigned int n = 0; n < m_nPlanes; n++) {
    steerFile << (counter * 3 + 1) << " " << m_xPosDet[n] << " 0.0" << endl;
    steerFile << (counter * 3 + 2) << " " << m_yPosDet[n] << " 0.0" << endl;
    steerFile << (counter * 3 + 3) << " " << m_gammaPosDet[n] << " 0.0" << endl;
    counter++;

  } // end loop over all planes

  steerFile << endl;
  steerFile << endl;
  steerFile << "method inversion 10 0.001" << endl;
  steerFile << endl;
  steerFile << "histprint" << endl;
  steerFile << endl;
  steerFile << "end" << endl;

  steerFile.close();

  streamlog_out ( MESSAGE5 ) << "File " << path << " written." << endl;
}
