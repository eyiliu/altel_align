#include "getopt.h"
#include "mysystem.hh"
#include "myrapidjson.h"

#include <iostream>
#include <algorithm>
#include <regex>

#include <Eigen/Geometry>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH2.h>
#include <TImage.h>
#include <TCanvas.h>


void FitTrack(unsigned int nMeasures,
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
                          std::vector<double>& xResidMeasure,
                          std::vector<double>& yResidMeasure);


static const std::string help_usage = R"(
Usage:
  -help                    help message
  -verbose                 verbose flag
  -hitFile         [PATH]   path to hit file  (input)
  -geometryFile    [PATH]   geometry input file (input)
  -rootFile        [PATH]   path to root file (output)
  -maxEventNumber   [int]    max number of events to be processed
  -maxTrackNumber   [int]    max number of tracks to be processed

examples:
../bin/LineFit_main -hitFile /work/data/TB2006/alpide_200629033515.json  -geo align_313_geo.json -r linefit.root -maxT 1000
)";

int main(int argc, char *argv[]) {
  int do_help = false;
  int do_verbose = false;
  struct option longopts[] = {{"help", no_argument, &do_help, 1},
                              {"verbose", no_argument, &do_verbose, 1},
                              {"geometryFile", required_argument, NULL, 'g'},
                              {"hitFile", required_argument, NULL, 'o'},
                              {"rootFile", required_argument, NULL, 'u'},
                              {"maxEventNumber", required_argument, NULL, 'm'},
                              {"maxTrackNumber", required_argument, NULL, 'n'},
                             {0, 0, 0, 0}};

  std::string geometryFile_path;
  std::string hitFile_path;
  std::string rootFile_path;
  size_t maxTrackNumber = -1;
  size_t maxEventNumber = -1;

  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL)) != -1) {
    switch (c) {
    case 'g':
      geometryFile_path = optarg;
      break;
    case 'o':
      hitFile_path = optarg;
      break;
    case 'u':
      rootFile_path = optarg;
      break;
    case 'm':
      maxEventNumber = std::stoull(optarg);
      break;
    case 'n':
      maxTrackNumber = std::stoull(optarg);
      break;
      /////generic part below///////////
    case 0: /* getopt_long() set a variable, just keep going */
      break;
    case 1:
      fprintf(stderr, "case 1\n");
      exit(1);
      break;
    case ':':
      fprintf(stderr, "case :\n");
      exit(1);
      break;
    case '?':
      fprintf(stderr, "case ?\n");
      exit(1);
      break;
    default:
      fprintf(stderr, "case default, missing branch in switch-case\n");
      exit(1);
      break;
    }
  }
  if(do_help){
    std::fprintf(stdout, "help message:\n");
    std::fprintf(stdout, "%s\n", help_usage.c_str());
    std::exit(1);
  }

  if (geometryFile_path.empty() ||
      hitFile_path.empty()||
      rootFile_path.empty()
    ) {
    std::fprintf(stderr, "%s\n", help_usage.c_str());
    std::exit(0);
  }

  std::fprintf(stdout, "\n");
  std::fprintf(stdout, "geometryFile: %s\n", geometryFile_path.c_str());
  std::fprintf(stdout, "hitFile:      %s\n", hitFile_path.c_str());
  std::fprintf(stdout, "rootFile:     %s\n", rootFile_path.c_str());
  std::fprintf(stdout, "\n");

  std::printf("--------read geo-----\n");
  std::string str_geo = JsonUtils::readFile(geometryFile_path.c_str());
  JsonDocument jsd_geo = JsonUtils::createJsonDocument(str_geo);
  if(jsd_geo.IsNull()){
    std::fprintf(stderr, "Geometry file <%s> does not contain any json objects.\n", geometryFile_path.c_str() );
    throw;
  }
  if( !jsd_geo.HasMember("geometry") ||
      !jsd_geo["geometry"].HasMember("detectors") ||
      !jsd_geo["geometry"]["detectors"].IsArray()){
    std::fprintf(stderr, "Geometry file <%s>: unknown json structure.\n", geometryFile_path.c_str() );
    JsonUtils::printJsonValue(jsd_geo, true);
    throw;
  }
  const auto &js_geo = jsd_geo["geometry"];
  const auto &js_dets = js_geo["detectors"];
  for(const auto& js_det : js_dets.GetArray()){
    JsonUtils::printJsonValue(js_det, false);
  }

  

  std::map<size_t, Eigen::Affine3d> trafo_dets;
  for(const auto& js_det : js_dets.GetArray()){
    size_t id = js_det["id"].GetUint();
    double cx = js_det["center"]["x"].GetDouble();
    double cy = js_det["center"]["y"].GetDouble();
    double cz = js_det["center"]["z"].GetDouble();
    double rx = js_det["rotation"]["x"].GetDouble();
    double ry = js_det["rotation"]["y"].GetDouble();
    double rz = js_det["rotation"]["z"].GetDouble();

    Eigen::AngleAxisd rotZ(rz, Eigen::Vector3d::UnitZ());
    Eigen::AngleAxisd rotY(ry, Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd rotX(rx, Eigen::Vector3d::UnitX());
    Eigen::Affine3d trafoGlobalFromMeas(Eigen::Affine3d::Identity());
    trafoGlobalFromMeas.rotate(rotX);
    trafoGlobalFromMeas.rotate(rotY);
    trafoGlobalFromMeas.rotate(rotZ);
    trafoGlobalFromMeas.translation() =  Eigen::Vector3d(cx, cy, cz);
    trafo_dets[id]=trafoGlobalFromMeas;
  }

  TTree tree("linefit", "linefit");
  size_t n_datapack_select_opt = 20000;
  JsonFileDeserializer jsf(hitFile_path);

  std::vector<size_t> idMeas;
  std::vector<double> xPosMeas;
  std::vector<double> yPosMeas;
  std::vector<double> zPosMeas;
  std::vector<double> xResolMeas;
  std::vector<double> yResolMeas;

  std::vector<double> xResidFit;
  std::vector<double> yResidFit;
  std::vector<double> xPosFit;
  std::vector<double> yPosFit;

  std::vector<size_t>* p_idMeas = &idMeas;
  std::vector<double>* p_xPosMeas = &xPosMeas;
  std::vector<double>* p_yPosMeas = &yPosMeas;
  std::vector<double>* p_zPosMeas = &zPosMeas;
  std::vector<double>* p_xResolMeas = &xResolMeas;
  std::vector<double>* p_yResolMeas = &yResolMeas;
  std::vector<double>* p_xResidFit = &xResidFit;
  std::vector<double>* p_yResidFit = &yResidFit;
  std::vector<double>* p_xPosFit = &xPosFit;
  std::vector<double>* p_yPosFit = &yPosFit;

  auto brId = tree.Branch("id", &p_idMeas);
  auto brXMeas = tree.Branch("xmeas", &p_xPosMeas);
  auto brYMeas = tree.Branch("ymeas", &p_yPosMeas);
  auto brXResid = tree.Branch("xresid", &p_xResidFit);
  auto brYResid = tree.Branch("yresid", &p_yResidFit);
  auto brXFit = tree.Branch("xfit", &p_xPosFit);
  auto brYFit = tree.Branch("yfit", &p_yPosFit);

  size_t nTracks = 0;
  size_t nEvents = 0;
  while (jsf) {
    if (nTracks >= maxTrackNumber || nEvents >= maxEventNumber ) {
      break;
    }
    auto evpack = jsf.getNextJsonDocument();
    if(evpack.IsNull()){
      std::fprintf(stdout, "reach null object, end of file\n");
      break;
    }
    const auto &layers = evpack["layers"];
    if (layers.Size() != trafo_dets.size()) {
      std::fprintf(stdout, "something wrong, skip envent %i\n", nEvents);
      nEvents++;
      continue;
    }

    idMeas.clear();
    xPosMeas.clear();
    yPosMeas.clear();
    zPosMeas.clear();
    xResolMeas.clear();
    yResolMeas.clear();

    xResidFit.clear();
    yResidFit.clear();
    xPosFit.clear();
    yPosFit.clear();

    for (const auto &layer : layers.GetArray()) {
      if(layer["hit"].Size() !=1){
        continue;
      }
      size_t id = layer["ext"].GetUint();
      for (const auto &hit : layer["hit"].GetArray()) {
        double x_hit = hit["pos"][0].GetDouble() - 0.02924 * 1024 / 2.0;
        double y_hit = hit["pos"][1].GetDouble() - 0.02688 * 512 / 2.0;
        Eigen::Vector3d measPos(x_hit, y_hit, 0);
        Eigen::Vector3d globalPos =  trafo_dets.at(id) * measPos;
        idMeas.push_back(id);
        xPosMeas.push_back(globalPos(0));
        yPosMeas.push_back(globalPos(1));
        zPosMeas.push_back(globalPos(2));
        xResolMeas.push_back(0.02);
        yResolMeas.push_back(0.02);
      }
    }

    if(idMeas.size() != layers.Size()){
      std::fprintf(stdout, "skip incompleted track in event %i\n", nEvents);
      nEvents++;
      continue;
    }
    //
    double xOriLine;
    double yOriLine;
    double xAngleLine;
    double yAngleLine;
    double xChisqLine;
    double yChisqLine;
    FitTrack(idMeas.size(), xPosMeas, yPosMeas, zPosMeas, xResolMeas, yResolMeas,
             xOriLine, yOriLine, xAngleLine, yAngleLine, xChisqLine, yChisqLine,
             xResidFit, yResidFit);
    nEvents++;
    nTracks++;

    for(int i=0; i<idMeas.size(); i++){
      xPosFit.push_back(xPosMeas[i]-xResidFit[i]);
      yPosFit.push_back(yPosMeas[i]-yResidFit[i]);
    }

    tree.Fill();
  }

  TFile tfile(rootFile_path.c_str(),"recreate");
  tree.Write();
  tfile.Close();
  return 0;
}




void FitTrack(unsigned int nMeasures,
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
              std::vector<double>& xResidMeasure,
              std::vector<double>& yResidMeasure){
  unsigned int nPlanesFit = nMeasures;
  const double * xPosFit = xPosMeasure.data();
  const double * yPosFit = yPosMeasure.data();
  const double * zPosFit = zPosMeasure.data();;
  const double * xResFit = xResolMeasure.data();
  const double * yResFit = yResolMeasure.data();

  xResidMeasure.resize(nPlanesFit);
  yResidMeasure.resize(nPlanesFit);
  double *residXFit = xResidMeasure.data();
  double *residYFit = yResidMeasure.data();

  // std::cout<< "------------"<<std::endl;
  // std::cout<< nMeasures<<std::endl;
  // std::cout<< xPosMeasure.size()<<std::endl;
  // std::cout<< yPosMeasure.size()<<std::endl;
  // std::cout<< zPosMeasure.size()<<std::endl;
  // std::cout<< xResolMeasure.size()<<std::endl;
  // std::cout<< yResolMeasure.size()<<std::endl;
  // std::cout<< xResidMeasure.size()<<std::endl;
  // std::cout<< yResidMeasure.size()<<std::endl;
  // std::cout<< "============="<<std::endl;


  float S1[2]   = {0,0};
  float Sx[2]   = {0,0};
  float Xbar[2] = {0,0};

  std::vector<float> Zbar_X_v(nPlanesFit);
  std::vector<float> Zbar_Y_v(nPlanesFit);

  float * Zbar_X = Zbar_X_v.data();
  float * Zbar_Y = Zbar_Y_v.data();

  for (unsigned int counter = 0; counter < nPlanesFit; counter++){
    Zbar_X[counter] = 0.;
    Zbar_Y[counter] = 0.;
  }

  float Sy[2]     = {0,0};
  float Ybar[2]   = {0,0};
  float Sxybar[2] = {0,0};
  float Sxxbar[2] = {0,0};
  float A2[2]     = {0,0};

  // define S1
  for(unsigned int  counter = 0; counter < nPlanesFit; counter++ ){
    S1[0] = S1[0] + 1/pow(xResFit[counter],2);
    S1[1] = S1[1] + 1/pow(yResFit[counter],2);
  }

  // define Sx
  for(unsigned int  counter = 0; counter < nPlanesFit; counter++ ){
    Sx[0] = Sx[0] + zPosFit[counter]/pow(xResFit[counter],2);
    Sx[1] = Sx[1] + zPosFit[counter]/pow(yResFit[counter],2);
  }

  // define Xbar
  Xbar[0]=Sx[0]/S1[0];
  Xbar[1]=Sx[1]/S1[1];

  // coordinate transformation !! -> bar
  for(unsigned int  counter = 0; counter < nPlanesFit; counter++ ){
    Zbar_X[counter] = zPosFit[counter]-Xbar[0];
    Zbar_Y[counter] = zPosFit[counter]-Xbar[1];
  }

  // define Sy
  for(unsigned int  counter = 0; counter < nPlanesFit; counter++ ){
    Sy[0] = Sy[0] + xPosFit[counter]/pow(xResFit[counter],2);
    Sy[1] = Sy[1] + yPosFit[counter]/pow(yResFit[counter],2);
  }

  // define Ybar
  Ybar[0]=Sy[0]/S1[0];
  Ybar[1]=Sy[1]/S1[1];

  // define Sxybar
  for(unsigned int  counter = 0; counter < nPlanesFit; counter++ ){
    Sxybar[0] = Sxybar[0] + Zbar_X[counter] * xPosFit[counter]/pow(xResFit[counter],2);
    Sxybar[1] = Sxybar[1] + Zbar_Y[counter] * yPosFit[counter]/pow(yResFit[counter],2);
  }

  // define Sxxbar
  for(unsigned int  counter = 0; counter < nPlanesFit; counter++ ){
    Sxxbar[0] = Sxxbar[0] + Zbar_X[counter] * Zbar_X[counter]/pow(xResFit[counter],2);
    Sxxbar[1] = Sxxbar[1] + Zbar_Y[counter] * Zbar_Y[counter]/pow(yResFit[counter],2);
  }

  // define A2
  A2[0]=Sxybar[0]/Sxxbar[0];
  A2[1]=Sxybar[1]/Sxxbar[1];

  // Calculate chi sqaured
  // Chi^2 for X and Y coordinate for hits in all planes
  for(unsigned int  counter = 0; counter < nPlanesFit; counter++ ){
    xChisqLine += pow(  -zPosFit[counter]*A2[0] +  xPosFit[counter] - Ybar[0] + Xbar[0]*A2[0]  ,2) /
      pow(  xResFit[counter]  ,2);

    xChisqLine += pow(-zPosFit[counter]*A2[1]
                      +yPosFit[counter]-Ybar[1]+Xbar[1]*A2[1],2)/pow(yResFit[counter],2);
  }

  for(unsigned int counter = 0; counter < nPlanesFit; counter++ ) {
    residXFit[counter] = xPosFit[counter] - (Ybar[0]-Xbar[0]*A2[0]+zPosFit[counter]*A2[0]); // sign reverse?
    residYFit[counter] = yPosFit[counter] - (Ybar[1]-Xbar[1]*A2[1]+zPosFit[counter]*A2[1]);
  }

  // define angle
  xAngleLine = atan(A2[0]);
  yAngleLine = atan(A2[1]);

  xOriginLine = Ybar[0] - Xbar[0]*A2[0];
  yOriginLine = Ybar[1] - Xbar[1]*A2[1];
}
