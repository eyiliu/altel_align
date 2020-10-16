// from: Igor Rubinskiy, DESY <mailto:igorrubinsky@gmail.com>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>


#include "streamlog/streamlog.h"
#include "pstream.h"
#include "Mille.h"
#include "EUTelMille.h"

using namespace std;

EUTelMille::EUTelMille (){

  _iEvt = 0;

  // Initialise Mille statistics
  _nMilleDataPoints = 0;
  _nMilleTracks = 0;

  _waferResidX = new double[_nPlanes];
  _waferResidY = new double[_nPlanes];
  _waferResidZ = new double[_nPlanes];

  _xFitPos = new double[_nPlanes];
  _yFitPos = new double[_nPlanes];

  _telescopeResolX = new double[_nPlanes];
  _telescopeResolY = new double[_nPlanes];
  _telescopeResolZ = new double[_nPlanes];

  _mille = new Mille(_binaryFilename.c_str());

}


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

  for (unsigned int help = 0; help < nPlanesFitter; help++) {


    xPosFit[nPlanesFit] = xPosFitter[help];
    yPosFit[nPlanesFit] = yPosFitter[help];
    zPosFit[nPlanesFit] = zPosFitter[help];
    xResFit[nPlanesFit] = xResFitter[help];
    yResFit[nPlanesFit] = yResFitter[help];
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


void EUTelMille::processEvent () {
  // fill resolution arrays
  for (size_t help = 0; help < _nPlanes; help++) {
    _telescopeResolX[help] = _telescopeResolution;
    _telescopeResolY[help] = _telescopeResolution;
  }

  int _nTracks = 0;

  std::vector<IntVec > indexarray;

  streamlog_out( DEBUG5 ) << "Event #" << _iEvt << std::endl;



  // TODO: parepare tracks
  // set _nTracks
  // i = trackN, j = planeN
  // _xPos[i][j] = _allHitsArray[j][indexarray[i][j]].measuredX;
  // _yPos[i][j] = _allHitsArray[j][indexarray[i][j]].measuredY;
  // _zPos[i][j] = _allHitsArray[j][indexarray[i][j]].measuredZ;

  streamlog_out( DEBUG5 ) << "Track finder found " << _nTracks << std::endl;


  // Perform fit for all found track candidates
  // ------------------------------------------
  double Chiquare[2] = {0,0};
  double angle[2] = {0,0};

  // loop over all track candidates
  for (int track = 0; track < _nTracks; track++) {
    _xPosHere = new double[_nPlanes];
    _yPosHere = new double[_nPlanes];
    _zPosHere = new double[_nPlanes];

    for (unsigned int help = 0; help < _nPlanes; help++) {
      _xPosHere[help] = _xPos[track][help];
      _yPosHere[help] = _yPos[track][help];
      _zPosHere[help] = _zPos[track][help];
    }

    Chiquare[0] = 0.0;
    Chiquare[1] = 0.0;

    streamlog_out ( MESSAGE1 ) << "Adding track using the following coordinates: ";

    // loop over all planes
    for (unsigned int help = 0; help < _nPlanes; help++) 
    {
      streamlog_out ( MESSAGE1 ) << 
        std::endl << " not Excluded @ " << help  << "["<<_nPlanes << "] " <<_xPosHere[help] << " " << _yPosHere[help] << " " << _zPosHere[help] ;
      streamlog_out ( MESSAGE1 ) << std::endl;

    } // end loop over all planes

    streamlog_out ( MESSAGE1 ) << endl;

    streamlog_out(MESSAGE1) << " AlignMode = " << _alignMode  << std::endl;
    // Calculate residuals
    FitTrack(_nPlanes,
             _xPosHere,
             _yPosHere,
             _zPosHere,
             _telescopeResolX,
             _telescopeResolY,
             Chiquare,
             _waferResidX,
             _waferResidY,
             angle);


    // Add track to Millepede
    // ---------------------------

    // Easy case: consider only shifts
    if (_alignMode == 2) {

      const int nLC = 4; // number of local parameters
      const int nGL = (_nPlanes) * 2; // number of global parameters

      float sigma = _telescopeResolution;

      float *derLC = new float[nLC]; // array of derivatives for local parameters
      float *derGL = new float[nGL]; // array of derivatives for global parameters

      int *label = new int[nGL]; // array of labels

      float residual;

      // create labels
      for (int help = 0; help < nGL; help++) {
        label[help] = help + 1;
      }

      for (int help = 0; help < nGL; help++) {
        derGL[help] = 0;
      }

      for (int help = 0; help < nLC; help++) {
        derLC[help] = 0;
      }

      // loop over all planes
      for (unsigned int help = 0; help < _nPlanes; help++) {
        int helphelp = help; // index of plane after

        derGL[((helphelp * 2) + 0)] = -1;
        derLC[0] = 1;
        derLC[2] = _zPosHere[help];
        residual = _waferResidX[help];
        sigma    = _resolutionX[help];
        _mille->mille(nLC,derLC,nGL,derGL,label,residual,sigma);

        derGL[((helphelp * 2) + 0)] = 0;
        derLC[0] = 0;
        derLC[2] = 0;

        derGL[((helphelp * 2) + 1)] = -1;
        derLC[1] = 1;
        derLC[3] = _zPosHere[help];
        residual = _waferResidY[help];
        sigma    = _resolutionY[help];
        _mille->mille(nLC,derLC,nGL,derGL,label,residual,sigma);

        derGL[((helphelp * 2) + 1)] = 0;
        derLC[1] = 0;
        derLC[3] = 0;

        _nMilleDataPoints++;
      } // end loop over all planes

      delete [] derLC;
      delete [] derGL;
      delete [] label;

      // Slightly more complicated: add rotation around the z axis
    } else if (_alignMode == 1) {

      const int nLC = 4; // number of local parameters
      const int nGL = _nPlanes * 3; // number of global parameters

      float sigma = _telescopeResolution;

      float *derLC = new float[nLC]; // array of derivatives for local parameters
      float *derGL = new float[nGL]; // array of derivatives for global parameters

      int *label = new int[nGL]; // array of labels

      float residual;

      // create labels
      for (int help = 0; help < nGL; help++) {
        label[help] = help + 1;
      }

      for (int help = 0; help < nGL; help++) {
        derGL[help] = 0;
      }

      for (int help = 0; help < nLC; help++) {
        derLC[help] = 0;
      }

      // loop over all planes
      for (unsigned int help = 0; help < _nPlanes; help++) {

        int helphelp = help;

        derGL[((helphelp * 3) + 0)] = -1;
        derGL[((helphelp * 3) + 2)] = _yPosHere[help];
        derLC[0] = 1;
        derLC[2] = _zPosHere[help];
        residual = _waferResidX[help];
        sigma    = _resolutionX[help];
        _mille->mille(nLC,derLC,nGL,derGL,label,residual,sigma);

        derGL[((helphelp * 3) + 0)] = 0;
        derGL[((helphelp * 3) + 2)] = 0;
        derLC[0] = 0;
        derLC[2] = 0;

        derGL[((helphelp * 3) + 1)] = -1;
        derGL[((helphelp * 3) + 2)] = -1 * _xPosHere[help];
        derLC[1] = 1;
        derLC[3] = _zPosHere[help];
        residual = _waferResidY[help];
        sigma    = _resolutionY[help];
        _mille->mille(nLC,derLC,nGL,derGL,label,residual,sigma);

        derGL[((helphelp * 3) + 1)] = 0;
        derGL[((helphelp * 3) + 2)] = 0;
        derLC[1] = 0;
        derLC[3] = 0;

        _nMilleDataPoints++;


      } // end loop over all planes

      // clean up

      delete [] derLC;
      delete [] derGL;
      delete [] label;

    }
    else if (_alignMode == 3) {
      const int nLC = 4; // number of local parameters
      const int nGL = _nPlanes * 6; // number of global parameters

      float *derLC = new float[nLC]; // array of derivatives for local parameters
      float *derGL = new float[nGL]; // array of derivatives for global parameters

      int *label = new int[nGL]; // array of labels

      float residual;

      // create labels
      for (int help = 0; help < nGL; help++)
      {
        label[help] = help + 1;
      }

      for (int help = 0; help < nGL; help++)
      {
        derGL[help] = 0;
      }

      for (int help = 0; help < nLC; help++)
      {
        derLC[help] = 0;
      }

      // loop over all planes
      for (unsigned int help = 0; help < _nPlanes; help++)
      {
        double sigmax = _resolutionX[help];
        double sigmay = _resolutionY[help];
        double sigmaz = _resolutionZ[help];

        int helphelp = help; // index of plane after

        for (int i = 0; i < nGL; i++ )
        {
          derGL[i] = 0.000;
        }

        for (int i = 0; i < nLC; i++ ) 
        {
          derLC[i] = 0.000;
        }

        double x_sensor = 0.;
        double y_sensor = 0.;
        double z_sensor = 0.;

        // int sensorID = _sensorIDVec[help];

        // //TODO:!!!!!  YI: store geo as data member , and retrived there
        // x_sensor = geo::gGeometry().siPlaneXPosition(sensorID); 
        // y_sensor = geo::gGeometry().siPlaneYPosition(sensorID); 
        // z_sensor = geo::gGeometry().siPlaneZPosition(sensorID);
        // //


// track model : fit-reco => 
//   (a_X*x+b_X, a_Y*y+b_Y)  :   /  1   -g    b \   / x          \   :: shouldn't it be x-xcenter-of-the-sensor ??
//                           : - |  g    1   -a |   | y          |   ::  and y-ycenter-of-the-sensor ??   
//                           :   \ -b    a    1 /   \ z-z_sensor /   ::  == z-zcenter-of-the-sensor  ?? (already)
// make angles sings consistent with X->Y->Z->X rotations.
// correct likewise all matrices in ApplyAlignment processor
// Igor Rubinsky 09-10-2011
//
        // shift in X
        derGL[((helphelp * 6) + 0)] = -1.0;                                 // dx
        derGL[((helphelp * 6) + 1)] =  0.0;                                 // dy
        derGL[((helphelp * 6) + 2)] =  0.0;                                 // dz
        // rotation in ZY ( alignment->getAlpfa() )
        derGL[((helphelp * 6) + 3)] =      0.0;                             // alfa  - ZY :: Y->Z
        derGL[((helphelp * 6) + 4)] = -1.0*(_zPosHere[help] - z_sensor);    // beta  - ZX :: Z->X
        derGL[((helphelp * 6) + 5)] =  1.0*(_yPosHere[help] - y_sensor);                 // gamma - XY :: X->Y

        derLC[0] = 1.0;
        derLC[1] = 0.0;
        derLC[2] = _zPosHere[help] + _waferResidZ[help];
        derLC[3] = 0.0;

        residual = _waferResidX[help];
        _mille->mille(nLC,derLC,nGL,derGL,label,residual,sigmax);


        // shift in Y
        derGL[((helphelp * 6) + 0)] =  0.0;
        derGL[((helphelp * 6) + 1)] = -1.0;
        derGL[((helphelp * 6) + 2)] =  0.0;
        // rotation in ZX
        derGL[((helphelp * 6) + 3)] =  1.0*(_zPosHere[help] - z_sensor);
        derGL[((helphelp * 6) + 4)] =      0.0            ;
        derGL[((helphelp * 6) + 5)] = -1.0*(_xPosHere[help] - x_sensor);

        derLC[0] = 0.0;
        derLC[1] = 1.0;
        derLC[2] = 0.0;
        derLC[3] = _zPosHere[help] + _waferResidZ[help];

        residual = _waferResidY[help];
        _mille->mille(nLC,derLC,nGL,derGL,label,residual,sigmay);
        // shift in Z
        derGL[((helphelp * 6) + 0)] =  0.0;
        derGL[((helphelp * 6) + 1)] =  0.0;
        derGL[((helphelp * 6) + 2)] = -1.0;
        // rotation in XY
        derGL[((helphelp * 6) + 3)] = -1.0*(_yPosHere[help]-y_sensor);
        derGL[((helphelp * 6) + 4)] =  1.0*(_xPosHere[help]-x_sensor);
        derGL[((helphelp * 6) + 5)] =  0.0;

        derLC[0] = 0.0;
        derLC[1] = 0.0;
        derLC[2] = _xPosHere[help] + _waferResidX[help];
        derLC[3] = _yPosHere[help] + _waferResidY[help];

        residual = _waferResidZ[help];
        _mille->mille(nLC,derLC,nGL,derGL,label,residual,sigmaz);
        _nMilleDataPoints++;

      } // end loop over all planes

      // clean up
      delete [] derLC;
      delete [] derGL;
      delete [] label;
    } else {
      streamlog_out ( ERROR2 ) << _alignMode << " is not a valid mode. Please choose 1,2 or 3." << endl;
    }

    // end local fit
    _mille->end();

    _nMilleTracks++;

    // Fill histograms for individual tracks
    // -------------------------------------

    // clean up
    delete [] _zPosHere;
    delete [] _yPosHere;
    delete [] _xPosHere;

  } // end loop over all track candidates

  streamlog_out ( MESSAGE1 ) << "Finished fitting tracks in event " << _iEvt << endl;
  // count events
  ++_iEvt;

}


void EUTelMille::end() {

  delete [] _telescopeResolY;
  delete [] _telescopeResolX;
  delete [] _telescopeResolZ;
  delete [] _yFitPos;
  delete [] _xFitPos;
  delete [] _waferResidY;
  delete [] _waferResidX;
  delete [] _waferResidZ;

  // close the output file
  delete _mille;

  streamlog_out ( MESSAGE4 ) << endl << "Generating the steering file for the pede program..." << endl;

  double *meanX = new double[_nPlanes];
  double *meanY = new double[_nPlanes];
  double *meanZ = new double[_nPlanes];

  // loop over all detector planes
  for(unsigned int iDetector = 0; iDetector < _nPlanes; iDetector++ ) {

    // TODO: YI 
    // meanX[iDetector] = residx_histo->mean(); // TODO, hist is removed
    // meanY[iDetector] = residy_histo->mean();
    // meanZ[iDetector] = residz_histo->mean();

    
  } // end loop over all detector planes

  ofstream steerFile;
  steerFile.open(_pedeSteerfileName.c_str());

  if (steerFile.is_open()) {

    unsigned int firstnotexcl = _nPlanes;
    unsigned int lastnotexcl = 0;

    // loop over all planes
    for (unsigned int help = 0; help < _nPlanes; help++) 
    {

      if (firstnotexcl > help) 
      {
        firstnotexcl = help;
      }

      if (lastnotexcl < help) 
      {
        lastnotexcl = help;
      }
    } // end loop over all planes

      // calculate average
    double averageX = (meanX[firstnotexcl] + meanX[lastnotexcl]) / 2;
    double averageY = (meanY[firstnotexcl] + meanY[lastnotexcl]) / 2;
    double averageZ = (meanZ[firstnotexcl] + meanZ[lastnotexcl]) / 2;

    steerFile << "Cfiles" << endl;
    steerFile << _binaryFilename << endl;
    steerFile << endl;

    steerFile << "Parameter" << endl;

    int counter = 0;

    // loop over all planes
    for (unsigned int help = 0; help < _nPlanes; help++) {

      bool fixed = false;
      for(size_t i = 0;i< _FixedPlanes.size(); i++)
      {
        if(_FixedPlanes[i] == static_cast< int >(help))
          fixed = true;
      }
          
      if( fixed || (_FixedPlanes.empty() && (help == firstnotexcl || help == lastnotexcl) ) )
      {
        if (_alignMode == 1) {
          steerFile << (counter * 3 + 1) << " 0.0 -1.0" << endl;
          steerFile << (counter * 3 + 2) << " 0.0 -1.0" << endl;
          steerFile << (counter * 3 + 3) << " 0.0 -1.0" << endl;
        } else if (_alignMode == 2) {
          steerFile << (counter * 2 + 1) << " 0.0 -1.0" << endl;
          steerFile << (counter * 2 + 2) << " 0.0 -1.0" << endl;
        } else if (_alignMode == 3) {
          steerFile << (counter * 6 + 1) << " 0.0 -1.0" << endl;
          steerFile << (counter * 6 + 2) << " 0.0 -1.0" << endl;
          steerFile << (counter * 6 + 3) << " 0.0 -1.0" << endl;
          steerFile << (counter * 6 + 4) << " 0.0 -1.0" << endl;
          steerFile << (counter * 6 + 5) << " 0.0 -1.0" << endl;
          steerFile << (counter * 6 + 6) << " 0.0 -1.0" << endl;
        }
              
      } else {
            
        if (_alignMode == 1) {

          if (_usePedeUserStartValues == 0) {
            steerFile << (counter * 3 + 1) << " " << (averageX - meanX[help]) << " 0.0" << endl;
            steerFile << (counter * 3 + 2) << " " << (averageY - meanY[help]) << " 0.0" << endl;
            steerFile << (counter * 3 + 3) << " " << " 0.0 0.0" << endl;
          } else {
            steerFile << (counter * 3 + 1) << " " << _pedeUserStartValuesX[help] << " 0.0" << endl;
            steerFile << (counter * 3 + 2) << " " << _pedeUserStartValuesY[help] << " 0.0" << endl;
            steerFile << (counter * 3 + 3) << " " << _pedeUserStartValuesGamma[help] << " 0.0" << endl;
          }

        } else if (_alignMode == 2) {

          if (_usePedeUserStartValues == 0) {
            steerFile << (counter * 2 + 1) << " " << (averageX - meanX[help]) << " 0.0" << endl;
            steerFile << (counter * 2 + 2) << " " << (averageY - meanY[help]) << " 0.0" << endl;
          } else {
            steerFile << (counter * 2 + 1) << " " << _pedeUserStartValuesX[help] << " 0.0" << endl;
            steerFile << (counter * 2 + 2) << " " << _pedeUserStartValuesY[help] << " 0.0" << endl;
          }

        } else if (_alignMode == 3) {
          if (_usePedeUserStartValues == 0)
          {
            if(_FixParameter[help] & (1 << 0))
              steerFile << (counter * 6 + 1) << " 0.0 -1.0" << endl;
            else
              steerFile << (counter * 6 + 1) << " " << (averageX - meanX[help]) << " 0.0" << endl;
                  
            if(_FixParameter[help] & (1 << 1))
              steerFile << (counter * 6 + 2) << " 0.0 -1.0" << endl;
            else
              steerFile << (counter * 6 + 2) << " " << (averageY - meanY[help]) << " 0.0" << endl;
                  
            if(_FixParameter[help] & (1 << 2))
              steerFile << (counter * 6 + 3) << " 0.0 -1.0" << endl;
            else
              steerFile << (counter * 6 + 3) << " " << (averageZ - meanZ[help]) << " 0.0" << endl;
                  
            if(_FixParameter[help] & (1 << 3))
              steerFile << (counter * 6 + 4) << " 0.0 -1.0" << endl;
            else
              steerFile << (counter * 6 + 4) << " 0.0 0.0" << endl;
                  
            if(_FixParameter[help] & (1 << 4))
              steerFile << (counter * 6 + 5) << " 0.0 -1.0" << endl;
            else
              steerFile << (counter * 6 + 5) << " 0.0 0.0" << endl;

            if(_FixParameter[help] & (1 << 5))
              steerFile << (counter * 6 + 6) << " 0.0 -1.0" << endl;
            else
              steerFile << (counter * 6 + 6) << " 0.0 0.0" << endl;
          }
          else
          {
            if(_FixParameter[help] & (1 << 0))
              steerFile << (counter * 6 + 1) << " 0.0 -1.0" << endl;
            else
              steerFile << (counter * 6 + 1) << " " << _pedeUserStartValuesX[help] << " 0.0" << endl;
                  
            if(_FixParameter[help] & (1 << 1))
              steerFile << (counter * 6 + 2) << " 0.0 -1.0" << endl;
            else
              steerFile << (counter * 6 + 2) << " " << _pedeUserStartValuesY[help] << " 0.0" << endl;
                  
            if(_FixParameter[help] & (1 << 2))
              steerFile << (counter * 6 + 3) << " 0.0 -1.0" << endl;
            else
              steerFile << (counter * 6 + 3) << " " << _pedeUserStartValuesZ[help] << " 0.0" << endl;
                  
            if(_FixParameter[help] & (1 << 3))
              steerFile << (counter * 6 + 4) << " 0.0 -1.0" << endl;
            else
              steerFile << (counter * 6 + 4) << " " << _pedeUserStartValuesAlpha[help] << " 0.0" << endl;
                  
            if(_FixParameter[help] & (1 << 4))
              steerFile << (counter * 6 + 5) << " 0.0 -1.0" << endl;
            else
              steerFile << (counter * 6 + 5) << " " << _pedeUserStartValuesBeta[help] << " 0.0" << endl;
                  
            if(_FixParameter[help] & (1 << 5))
              steerFile << (counter * 6 + 6) << " 0.0 -1.0" << endl;
            else
              steerFile << (counter * 6 + 6) << " " << _pedeUserStartValuesGamma[help] << " 0.0" << endl;
          }
        }
      }
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

    streamlog_out ( MESSAGE5 ) << "File " << _pedeSteerfileName << " written." << endl;

  } else {

    streamlog_out ( ERROR2 ) << "Could not open steering file." << endl;

  }
  //cleaning up
  delete [] meanX;
  delete [] meanY;
  delete [] meanZ;

  streamlog_out ( MESSAGE7 ) << "Number of data points used: " << _nMilleDataPoints << endl;
  streamlog_out ( MESSAGE7 ) << "Number of tracks used: " << _nMilleTracks << endl;

  // if running pede using the generated steering file

  // check if steering file exists

  std::string command = "pede " + _pedeSteerfileName;

  streamlog_out ( MESSAGE5 ) << "Starting pede...: " << command.c_str() << endl;

  bool encounteredError = false;
  // run pede and create a streambuf that reads its stdout and stderr
  redi::ipstream pede( command.c_str(), redi::pstreams::pstdout|redi::pstreams::pstderr );

  if (!pede.is_open()) {
    streamlog_out( ERROR5 ) << "Pede cannot be executed: command not found in the path" << endl;
    encounteredError = true;
  } else {
    // output multiplexing: parse pede output in both stdout and stderr and echo messages accordingly
    char buf[1024];
    std::streamsize n;
    std::stringstream pedeoutput; // store stdout to parse later
    std::stringstream pedeerrors;
    bool finished[2] = { false, false };
    while (!finished[0] || !finished[1])
    {
      if (!finished[0])
      {
        while ((n = pede.err().readsome(buf, sizeof(buf))) > 0){
          streamlog_out( ERROR5 ).write(buf, n).flush();
          string error (buf, n);
          pedeerrors << error;
          encounteredError = true;
        }
        if (pede.eof())
        {
          finished[0] = true;
          if (!finished[1])
            pede.clear();
        }
      }

      if (!finished[1])
      {
        while ((n = pede.out().readsome(buf, sizeof(buf))) > 0){
          streamlog_out( MESSAGE4 ).write(buf, n).flush();
          string output (buf, n);
          pedeoutput << output;
        }
        if (pede.eof())
        {
          finished[1] = true;
          if (!finished[0])
            pede.clear();
        }
      }
    }

    // pede does not return exit codes on some errors (in V03-04-00)
    // check for some of those here by parsing the output
    {
      const char * pch = strstr(pedeoutput.str().data(),"Too many rejects");
      if (pch){
        streamlog_out ( ERROR5 ) << "Pede stopped due to the large number of rejects. " << endl;
        encounteredError = true;
      }
    }
	
    {
      const char* pch0 = strstr(pedeoutput.str().data(),"Sum(Chi^2)/Sum(Ndf) = ");
      if (pch0 != 0){
        streamlog_out ( DEBUG5 ) << " Parsing pede output for final chi2/ndf result.. " << endl;
        // search for the equal sign after which the result for chi2/ndf is stated within the next 80 chars 
        // (with offset of 22 chars since pch points to beginning of "Sum(..." string just found)
        char* pch = (char*)((memchr (pch0+22, '=', 180)));
        if (pch!=NULL){
          char str[16];
          // now copy the numbers after the equal sign
          strncpy ( str, pch+1, 15 );
          str[15] = '\0';   /* null character manually added */
          // monitor the chi2/ndf in CDash when running tests
          streamlog_out ( MESSAGE6 ) << "Final Sum(Chi^2)/Sum(Ndf) = " << str << endl;
        }	    
      }
    }

    // wait for the pede execution to finish
    pede.close();

    // reading back the millepede.res file and getting the
    // results.
    string millepedeResFileName = "millepede.res";

    streamlog_out ( MESSAGE6 ) << "Reading back the " << millepedeResFileName << endl;

    // open the millepede ASCII output file
    ifstream millepede( millepedeResFileName.c_str() );


    /*
      if ( millepede.bad() || !millepede.is_open() )
      {
      streamlog_out ( ERROR4 ) << "Error opening the " << millepedeResFileName << endl
      << "The alignment slcio file cannot be saved" << endl;
      }
      else 
      {
      vector<double > tokens;
      stringstream tokenizer;
      string line;

      // get the first line and throw it away since it is a
      // comment!
      getline( millepede, line );

      int counter = 0;

      while ( ! millepede.eof() ) {

      EUTelAlignmentConstant * constant = new EUTelAlignmentConstant;

      bool goodLine = true;
      unsigned int numpars = 0;
      if(_alignMode != 3)
      numpars = 3;
      else
      numpars = 6;

      for ( unsigned int iParam = 0 ; iParam < numpars ; ++iParam ) 
      {
      getline( millepede, line );

      if ( line.empty() ) {
      goodLine = false;
      continue;
      }

      tokens.clear();
      tokenizer.clear();
      tokenizer.str( line );

      double buffer;
      // // check that all parts of the line are non zero
      while ( tokenizer >> buffer ) {
      tokens.push_back( buffer ) ;
      }

      if ( ( tokens.size() == 3 ) || ( tokens.size() == 6 ) || (tokens.size() == 5) ) {
      goodLine = true;
      } else goodLine = false;

      bool isFixed = ( tokens.size() == 3 );
      if(_alignMode != 3)
      {
      if ( iParam == 0 ) {
      constant->setXOffset( tokens[1] / 1000. );
      if ( ! isFixed ) constant->setXOffsetError( tokens[4] / 1000. ) ;
      }
      if ( iParam == 1 ) {
      constant->setYOffset( tokens[1] / 1000. ) ;
      if ( ! isFixed ) constant->setYOffsetError( tokens[4] / 1000. ) ;
      }
      if ( iParam == 2 ) {
      constant->setGamma( tokens[1]  ) ;
      if ( ! isFixed ) constant->setGammaError( tokens[4] ) ;
      }
      }
      else
      {
      if ( iParam == 0 ) {
      constant->setXOffset( tokens[1] / 1000. );
      if ( ! isFixed ) constant->setXOffsetError( tokens[4] / 1000. ) ;                    
      }
      if ( iParam == 1 ) {
      constant->setYOffset( tokens[1] / 1000. ) ;
      if ( ! isFixed ) constant->setYOffsetError( tokens[4] / 1000. ) ;
      }
      if ( iParam == 2 ) {
      constant->setZOffset( tokens[1] / 1000. ) ;
      if ( ! isFixed ) constant->setZOffsetError( tokens[4] / 1000. ) ;
      }
      if ( iParam == 3 ) {
      constant->setAlpha( tokens[1]  ) ;
      if ( ! isFixed ) constant->setAlphaError( tokens[4] ) ;
      }
      if ( iParam == 4 ) {
      constant->setBeta( tokens[1]  ) ;
      if ( ! isFixed ) constant->setBetaError( tokens[4] ) ;
      }
      if ( iParam == 5 ) {
      constant->setGamma( tokens[1]  ) ;
      if ( ! isFixed ) constant->setGammaError( tokens[4] ) ;
      }

      }
      }

      // right place to add the constant to the collection
      if ( goodLine  ) {
      constant->setSensorID( _orderedSensorID.at( counter ) );
      ++ counter;
      streamlog_out ( MESSAGE0 ) << (*constant) << endl;
      }
      else delete constant;
      }
      }
    */
    millepede.close();
  }

  streamlog_out ( MESSAGE2 ) << endl;
  streamlog_out ( MESSAGE2 ) << "Successfully finished" << endl;
}
