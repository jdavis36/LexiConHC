#include <cassert>
#include <exception>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "LexiConHCCouplings.h"
#include "LexiConHCTranslator.h"
#include "LexiConHCHelperFunctions.h"


using namespace std;
using namespace LexiConHCHelperFunctions;
using namespace LexiConHCIOHelpers;
using namespace LexiConHCCouplings;


LexiConHCTranslator::LexiConHCTranslator(LexiConHCOptionParser const& opts_) :
  opts(opts_)
{
  translate();
}
void LexiConHCTranslator::translate(){
  auto const& input_flags = opts.getInputFlags();
  auto const& input_parameters = opts.getInputParameters();
  auto const& input_couplings = opts.getInputCouplings();
  auto const& basis_input = opts.getInputBasis();
  auto const& basis_output = opts.getOutputBasis();

  // Get the translation matrix
  std::vector<std::vector<double>> tmatrix = getTranslationMatrix(basis_input, basis_output, input_parameters);
  // Get the input vector
  std::vector< std::pair<double, double> > vinput = getOrderedInputCouplings(basis_input, input_flags, input_parameters, input_couplings);
  // Assign the output vector
  std::vector< std::pair<double, double> > voutput(tmatrix.size(), std::pair<double, double>(0, 0));
  for (int i = 0; i < tmatrix.size(); i++)
      {
        for (int j = 0; j < tmatrix[i].size(); j++)
        {
            cout << tmatrix[i][j] << ",";
        }
        cout << "\n";
      }
  for (int i = 0; i < vinput.size(); i++)
  {
	cout << "(" << vinput[i].first << "," << vinput[i].second << ")";
  }
  cout <<"\n";
  for (int i = 0; i < voutput.size(); i++)
  {
        cout << "(" << voutput[i].first << "," << voutput[i].second << ")";
  }
  // Do the matrix multiplication
  for (size_t i=0; i<tmatrix.size(); i++){
    auto const& trow = tmatrix.at(i);
    for (size_t j=0; j<trow.size(); j++){
      voutput.at(i).first += trow.at(j) * vinput.at(j).first;
      voutput.at(i).second += trow.at(j) * vinput.at(j).second;
    }
  }
  // Fix output by subracting by constant vector 
  if (basis_input == bAmplitude_JHUGen && basis_output == bEFT_HiggsBasis){
     voutput.at(0).first -= 1;
     voutput.at(4).first -= 1;
  }
  else if (basis_input == bEFT_JHUGen && basis_output == bHiggsBasis){
     voutput.at(0).first -= 1;
     voutput.at(4).first -= 1;
  }
  else if (basis_input == bEFT_JHUGen && basis_output == bEFT_HiggsBasis){
     voutput.at(0).first -= 1;
  }
  else if (basis_input == bHiggsBasis && basis_output == bAmplitude_JHUGen){
     voutput.at(0).first += 2;
     voutput.at(4).first += 2;
  }
  else if (basis_input == bEFT_HiggsBasis && basis_output == bAmplitude_JHUGen){
     voutput.at(0).first += 2;
     voutput.at(4).first += 2;
  }
  // Set the results
  interpretOutputCouplings(basis_output, input_flags, input_parameters, voutput);
}

std::vector<std::vector<double>> LexiConHCTranslator::getTranslationMatrix(
  LexiConHCIOHelpers::IOBasisType const& basis_input, LexiConHCIOHelpers::IOBasisType const& basis_output,
  std::unordered_map<std::string, double> const& input_parameters
) const{
  double e; getValueWithDefault<std::string, double>(input_parameters, "e", e, DEFVAL_e);
  double gs; getValueWithDefault<std::string, double>(input_parameters, "gs", gs, DEFVAL_gs);
  double sw; getValueWithDefault<std::string, double>(input_parameters, "sin2ThetaW", sw, DEFVAL_SW);
  double cw; getValueWithDefault<std::string, double>(input_parameters, "cos2ThetaW", cw, DEFVAL_CW);
  double MZ; getValueWithDefault<std::string, double>(input_parameters, "MZ", MZ, DEFVAL_MZ);
  double MW; getValueWithDefault<std::string, double>(input_parameters, "MW", MW, DEFVAL_MW);
  double Lambda_z1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_z1", Lambda_z1, DEFVAL_LAMBDA_VI);
  double Lambda_w1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_w1", Lambda_w1, DEFVAL_LAMBDA_VI);
  double Lambda_zgs1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_zgs1", Lambda_zgs1, DEFVAL_LAMBDA_VI);
  std::vector<std::vector<double>> res;

  // This is where the translation matrixes should be coded carefully.
  // See also the getOrderedInputCouplings and interpretOutputCouplings functions for what the input/output vectors expect for JHUGen conventions (e.g. ghz1_prime2 scaled already by MZ^2/L1ZZ^2, so no need to put that for example).
  if (basis_input == bAmplitude_JHUGen){
    if (basis_output == bAmplitude_JHUGen){
      res.assign(nAmplitude_JHUGen_CouplingTypes, std::vector<double>(nAmplitude_JHUGen_CouplingTypes, 0));//Assign Correct Dimensions of translation matrix
      for (size_t i=0; i<(size_t) nAmplitude_JHUGen_CouplingTypes; i++) res.at(i).at(i)=1;
    }
    else if (basis_output == bEFT_JHUGen){
      res.assign(nEFT_JHUGen_CouplingTypes, std::vector<double>(nAmplitude_JHUGen_CouplingTypes, 0));
      
    }
    else if (basis_output == bHiggsBasis){
      res.assign(nEFT_HiggsBasis_CouplingTypes, std::vector<double>(nAmplitude_JHUGen_CouplingTypes, 0));
      res = std::vector<std::vector<double>> { 
	{0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,(-2*pow(cw,2)*pow(sw,2))/pow(e,2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
   	{0,(pow(MZ,2)*pow(sw,2))/(pow(e,2)*pow(Lambda_z1,2)),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
   	{0,0,0,(-2*pow(cw,2)*pow(sw,2))/pow(e,2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
   	{0,0,0,0,0,0,(-2*pow(sw,2))/pow(e,2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,(-2*pow(cw,2)*pow(sw,2))/(pow(e,2)*pow(Lambda_w1,2)),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
   	{0,0,0,0,0,0,0,(-2*pow(sw,2))/pow(e,2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,(-2*(sw) * cw /pow(e,2)),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
   	{0,0,0,0,0,0,0,0,0,0,(-2*(sw) * cw /pow(e,2)),0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,(cw*pow(MZ,2)*sw)/(pow(e,2)*pow(Lambda_zgs1,2)),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
   	{0,0,0,0,0,0,0,0,0,0,0,0,-2/pow(e,2),0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,-2/pow(e,2),0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2/pow(gs,2),0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2/pow(gs,2),0,0,0,0,0,0,0,0,0,0} };
      
    }
    else if (basis_output == bEFT_HiggsBasis){
      res.assign(nEFT_HiggsBasis_CouplingTypes, std::vector<double>(nAmplitude_JHUGen_CouplingTypes, 0));
    /*  res = std::vector<std::vector<double>> {
        {0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,(-2*pow(DEFVAL_CW,2)*pow(DEFVAL_SW,2))/pow(DEFVAL_e,2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,(pow(DEFVAL_MZ,2)*pow(DEFVAL_SW,2))/(pow(DEFVAL_e,2)*pow(DEFVAL_LAMBDA_VI,2)),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,(-2*pow(DEFVAL_CW,2)*pow(DEFVAL_SW,2))/pow(DEFVAL_e,2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,(-2*pow(DEFVAL_SW,2))/pow(DEFVAL_e,2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,(-2*pow(DEFVAL_CW,2)*pow(DEFVAL_SW,2))/pow(DEFVAL_e,2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,(-2*pow(DEFVAL_SW,2))/pow(DEFVAL_e,2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,(-2*(DEFVAL_SW) * DEFVAL_CW /pow(DEFVAL_e,2)),0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,(-2*(DEFVAL_SW) * DEFVAL_CW /pow(DEFVAL_e,2)),0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,(DEFVAL_CW*pow(DEFVAL_MZ,2)*DEFVAL_SW)/(pow(DEFVAL_e,2)*pow(DEFVAL_LAMBDA_VI,2)),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,-2/pow(DEFVAL_e,2),0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,-2/pow(DEFVAL_e,2),0,0,0,0,0,0,0,0,0,0} }; */
    }
  }
  else if (basis_input == bEFT_JHUGen){
    if (basis_output == bAmplitude_JHUGen){
      res.assign(nAmplitude_JHUGen_CouplingTypes, std::vector<double>(nEFT_JHUGen_CouplingTypes, 0));
      res = std::vector<std::vector<double>> {
					     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                             {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                                             {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
                                             {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
                                             {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                             {0,1/pow(Lambda_z1,2) * pow(Lambda_w1,2)/(pow(cw,2) - pow(sw,2)),(-2*pow(sw,2))/pow(MZ,2) * pow(Lambda_w1,2)/(pow(cw,2) - pow(sw,2)),0,(2*sw)/cw *(pow(cw,2) - pow(sw,2))/pow(MZ,2)*pow(Lambda_w1,2)/(pow(cw,2) - pow(sw,2)),0,(2*pow(sw,2))/pow(MZ,2)*pow(Lambda_w1,2)/(pow(cw,2) - pow(sw,2)),0, 0, 0},
                                             {0, 0, pow(cw,2), 0, 2*sw*cw, 0, pow(sw,2), 0, 0, 0},
                                             {0, 0, 0, pow(cw,2), 0, 2*sw*cw, 0, pow(sw,2), 0, 0},
                                             {0, (2*sw*cw)/pow(Lambda_z1,2)*pow(Lambda_zgs1,2)/(pow(cw,2) - pow(sw,2)),(-2*sw*cw)/pow(MZ,2),pow(Lambda_zgs1,2)/(pow(cw,2) - pow(sw,2)),0,(2*(pow(cw,2) - pow(sw,2)))/pow(MZ,2) * pow(Lambda_zgs1,2)/(pow(cw,2) - pow(sw,2)),0,(2*sw*cw)/pow(MZ,2)*pow(Lambda_zgs1,2)/(pow(cw,2) - pow(sw,2)), 0, 0, 0},
                                             {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
                                             {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
                                             {0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                                             {0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
                                             {0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
                                             {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                                             {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                             {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                             {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                             {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                             {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                             {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                             {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		   			     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                             {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
   					     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                             {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    }
    else if (basis_output == bEFT_JHUGen){
      res.assign(nEFT_JHUGen_CouplingTypes, std::vector<double>(nEFT_JHUGen_CouplingTypes, 0));
      for (size_t i=0; i<(size_t) nEFT_JHUGen_CouplingTypes; i++) res.at(i).at(i)=1;
    }
    else if (basis_output == bHiggsBasis){
      res.assign(nHiggsBasis_CouplingTypes, std::vector<double>(nEFT_JHUGen_CouplingTypes, 0)); 
      res = std::vector<std::vector<double>>{
				      {1/2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
				      {0, 0, -((2*pow(cw,2)*pow(sw,2))/pow(e,2)), 0, 0, 0, 0, 0, 0 ,0},
				      {0, (pow(MZ,2)*pow(sw,2))/(pow(e,2)*pow(Lambda_z1,2)), 0, 0, 0, 0, 0, 0, 0, 0},
				      {0, 0, 0, -((2*pow(cw,2)*pow(sw,2))/pow(e,2)), 0, 0, 0, 0, 0, 0},
				      {1/2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
				      {0, 0, -(2*(pow(cw,2)*pow(sw,2))/pow(e,2)), 0, -((4*cw*pow(sw,3))/pow(e,2)), 0, -((2*pow(sw,4))/pow(e,2)), 0, 0, 0},
				      {0, (pow(MW,2)*pow(sw,2))/(pow(e,2)*Lambda_z1*(pow(cw,2) - pow(sw,2))), -((2*pow(MW,2)*pow(sw,4))/(pow(e,2)*pow(MZ,2)*(pow(cw,2) - pow(sw,2)))),0,(2*pow(MW,2)*pow(sw,3))/(cw*pow(e,2)*pow(MZ,2)), 0,(2*pow(MW,2)*pow(sw,4))/(pow(e,2)*pow(MZ,2)*(pow(cw,2) - pow(sw,2))), 0, 0, 0},
				      {0, 0, 0, -((2*pow(cw,2)*pow(sw,2))/pow(e,2)), 0,-((4*cw*pow(sw,3))/pow(e,2)), 0, -((2*pow(sw,4))/pow(e,2)), 0, 0},
				      {0, 0, 0, 0, -((2*cw*sw)/pow(e,2)), 0, 0, 0, 0, 0},
				      {0, 0, 0, 0, 0, -((2*cw*sw)/pow(e,2)), 0, 0, 0, 0},
				      {0, (2*pow(cw,2)*pow(MZ,2)*pow(sw,2))/(pow(e,2)*pow(Lambda_z1,2)*(pow(cw,2) - pow(sw,2))), -((2*pow(cw,2)*pow(sw,2))/(pow(e,2)*(pow(cw,2) - pow(sw,2)))),0,(2*cw*sw)/pow(e,2), 0, (2*pow(cw,2)*pow(sw,2))/(pow(e,2)*(pow(cw,2) - pow(sw,2))), 0, 0, 0},
				      {0, 0, 0, 0, 0, 0, -(2/pow(e,2)), 0, 0, 0},
				      {0, 0, 0, 0, 0, 0, 0, -(2/pow(e,2)), 0, 0},
				      {0, 0, 0, 0, 0, 0, 0, 0, -2/pow(gs,2), 0},
				      {0, 0, 0, 0, 0, 0, 0, 0, 0, -2/pow(gs,2)}};
      }
    else if (basis_output == bEFT_HiggsBasis){
      res.assign(nHiggsBasis_CouplingTypes, std::vector<double>(nEFT_JHUGen_CouplingTypes, 0));
      res = std::vector<std::vector<double>>{
				    	    {1/2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
					    {0, (pow(MZ,2)* pow(sw,2))/(pow(e,2)*pow(Lambda_z1,2)), 0, 0, 0, 0, 0, 0, 0, 0},
					    {0, 0, (-2*pow(sw,2)*pow(cw,2))/pow(e,2), 0, 0, 0, 0, 0, 0, 0},
					    {0, 0, 0, (-2*pow(sw,2)*pow(cw,2))/pow(e,2), 0, 0, 0, 0, 0, 0},
					    {0, 0, 0, 0, (-2*sw*cw)/pow(e,2), 0, 0, 0 , 0, 0},
					    {0, 0, 0, 0, 0, (-2*sw*cw)/pow(e,2), 0, 0, 0, 0},
					    {0, 0, 0, 0, 0, 0, (-2)/pow(e,2), 0, 0, 0},
					    {0, 0, 0, 0, 0, 0, 0, -2/pow(e,2), 0, 0},
					    {0, 0, 0, 0, 0, 0, 0, 0, -2/pow(gs,2), 0},
					    {0, 0, 0, 0, 0, 0, 0, 0, 0, -2/pow(gs,2)}};
  }
}
  else if (basis_input == bHiggsBasis){
    cout << "HERE" ;
    if (basis_output == bAmplitude_JHUGen){
      res.assign(nAmplitude_JHUGen_CouplingTypes, std::vector<double>(nHiggsBasis_CouplingTypes, 0));
        res = std::vector<std::vector<double>>{
				{2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  				{0, 0, (pow(e,2)*pow(Lambda_z1,2))/(pow(MZ,2)*pow(sw,2)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  				{0, -(pow(e,2)/(2*pow(cw,2) * pow(sw,2))), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  				{0, 0, 0, -(pow(e,2)/(2*pow(cw,2)*pow(sw,2))), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  				{0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  				{0, 0, 0, 0, 0, 0, (pow(e,2)*pow(Lambda_w1,2))/(pow(MW,2)*pow(sw,2)), 0, 0, 0, 0, 0, 0, 0, 0},
  				{0, 0, 0, 0, 0, -(pow(e,2)/(2*pow(sw,2))), 0, 0, 0, 0, 0, 0, 0, 0, 0},
  				{0, 0, 0, 0, 0, 0, 0, -(pow(e,2)/(2*pow(sw,2))), 0, 0, 0, 0, 0, 0, 0},
  				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (pow(e,2)*pow(Lambda_zgs1,2))/(cw*pow(MZ,2)*sw), 0, 0, 0, 0},
  				{0, 0, 0, 0, 0, 0, 0, 0, -(pow(e,2)/(2*cw*sw)), 0, 0, 0, 0, 0, 0},
  				{0, 0, 0, 0, 0, 0, 0, 0, 0, -(pow(e,2)/(2*cw*sw)), 0, 0, 0, 0, 0},
  				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(pow(e,2)/2), 0, 0, 0},
  				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(pow(e,2)/2), 0, 0},
  				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(pow(gs,2)/2), 0},
  				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(pow(gs,2)/2)},
  				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
 
    }
    else if (basis_output == bEFT_JHUGen){
      res.assign(nEFT_JHUGen_CouplingTypes, std::vector<double>(nHiggsBasis_CouplingTypes, 0));
      /*
 *       res = std::vector<std::vector<double>>{
 *               {  },
 *                       {  }
 *                             };
 *                                */
    }
    else if (basis_output == bEFT_HiggsBasis){
      res.assign(nEFT_HiggsBasis_CouplingTypes, std::vector<double>(nHiggsBasis_CouplingTypes, 0));
      
    }
    else if (basis_output == bHiggsBasis){
      res.assign(nHiggsBasis_CouplingTypes, std::vector<double>(nHiggsBasis_CouplingTypes, 0));
      for (size_t i=0; i<(size_t) nHiggsBasis_CouplingTypes; i++) res.at(i).at(i)=1;
    }
  }
  else if (basis_input == bEFT_HiggsBasis){
    if (basis_output == bAmplitude_JHUGen){
      res.assign(nAmplitude_JHUGen_CouplingTypes, std::vector<double>(nEFT_HiggsBasis_CouplingTypes, 0));
      res = std::vector<std::vector<double>>{
					{2, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					{0, (pow(e,2)*pow(Lambda_z1,2))/(pow(MZ,2)*pow(sw,2)), 0, 0, 0, 0, 0, 0, 0, 0}, 
					{0, 0, -(pow(e,2)/(2*pow(cw,2)*pow(sw,2))), 0, 0, 0, 0, 0, 0, 0},
					{0, 0, 0, -(pow(e,2))/(2*pow(cw,2)*pow(sw,2)), 0, 0, 0, 0, 0, 0}, 
					{2, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					{0, (pow(e,2)*pow(Lambda_w1,2))/(pow(MZ,2)*pow(sw,2)*(pow(cw,2) - pow(sw,2))),(pow(e,2)*pow(Lambda_w1,2))/(pow(cw,2)*pow(MZ,2)*(pow(cw,2) - pow(sw,2))), 0, -((pow(e,2)*pow(Lambda_w1,2))/(pow(cw,2)*pow(MZ,2))), 0, -((pow(e,2)*pow(Lambda_w1,2)*pow(sw,2))/(pow(MZ,2)*(pow(cw,2) - pow(sw,2)))), 0, 0, 0}, 
					{0, 0, -(pow(e,2)/(2*pow(sw,2))), 0, -pow(e,2), 0, -(1/2)*pow(e,2)*pow(sw,2), 0, 0, 0}, 
					{0, 0, 0, -(pow(e,2))/(2*pow(sw,2)), 0, -pow(e,2), 0, -(1/2)*pow(e,2)*pow(sw,2), 0, 0}, 
					{0, (2*cw*pow(e,2)*pow(Lambda_zgs1,2))/(pow(MZ,2)*sw*(pow(cw,2) - pow(sw,2))),(pow(e,2)*pow(Lambda_zgs1,2))/(cw*pow(MZ,2)*sw*(pow(cw,2) - pow(sw,2))), 0, -((pow(e,2)*pow(Lambda_zgs1,2))/(cw*pow(MZ,2)*sw)), 0, -((cw*pow(e,2)* pow(Lambda_zgs1,2)*sw)/(pow(MZ,2)*(pow(cw,2) - pow(sw,2)))), 0, 0, 0},
					{0, 0, 0, 0, -(pow(e,2)/(2*cw*sw)), 0, 0, 0, 0, 0}, 
					{0, 0, 0, 0, 0, -(pow(e,2)/(2*cw*sw)), 0, 0, 0, 0}, 
					{0, 0, 0, 0, 0, 0, -(pow(e,2)/2), 0, 0, 0},
					{0, 0, 0, 0, 0, 0, 0, -(pow(e,2)/2), 0, 0}, 
					{0, 0, 0, 0, 0, 0, 0, 0, -(pow(gs,2)/2), 0}, 
					{0, 0, 0, 0, 0, 0, 0, 0, 0, -(pow(gs,2)/2)}, 
					{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    }
    else if (basis_output == bEFT_JHUGen){
      res.assign(nEFT_JHUGen_CouplingTypes, std::vector<double>(nEFT_HiggsBasis_CouplingTypes, 0));
      res = std::vector<std::vector<double>>{
				{2, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
				{0, (pow(e,2)*pow(Lambda_z1,2))/(pow(MZ,2)*pow(sw,2)), 0, 0, 0, 0, 0, 0, 0, 0}, 
				{0, 0, -(pow(e,2))/(2*pow(cw,2)*pow(sw,2)), 0, 0, 0, 0, 0, 0, 0}, 
				{0, 0, 0, -pow(e,2)/(2*pow(cw,2)*pow(sw,2)), 0, 0, 0, 0, 0, 0}, 
				{0, 0, 0, 0, -(pow(e,2)/(2*cw*sw)), 0, 0, 0, 0, 0}, 
				{0, 0, 0, 0, 0, -(pow(e,2)/(2*cw*sw)), 0, 0, 0, 0}, 
				{0, 0, 0, 0, 0, 0, -(pow(e,2)/2), 0, 0, 0}, 
				{0, 0, 0, 0, 0, 0, 0, -(pow(e,2)/2), 0, 0}, 
				{0, 0, 0, 0, 0, 0, 0, 0, -(pow(gs,2)/2), 0}, 
				{0, 0, 0, 0, 0, 0, 0, 0, 0, -(pow(gs,2)/2)}};
    }
    else if (basis_output == bHiggsBasis){
      res.assign(nHiggsBasis_CouplingTypes, std::vector<double>(nEFT_HiggsBasis_CouplingTypes, 0));
      
         res = std::vector<std::vector<double>>{
					{1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					{0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, 
					{0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, 
					{0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, 
					{1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					{0, 0, 1, 0, 2*pow(sw,2), 0,pow(sw,4), 0, 0, 0}, 
					{0, pow(MW,2)/(pow(MZ,2)*(pow(cw,2) - pow(sw,2))),(pow(MW,2)*pow(sw,2))/(pow(cw,2)*pow(MZ,2)*(pow(cw,2) - pow(sw,2))), 0, -((pow(MW,2)*pow(sw,2))/(pow(cw,2)*pow(MZ,2))), 0, -((pow(MW,2)*pow(sw,4))/(pow(MZ,2)*(pow(cw,2) - pow(sw,2)))), 0, 0, 0}, 
					{0, 0, 0, 1, 0, 2*pow(sw,2),0, pow(sw,4), 0, 0}, 
					{0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, 
					{0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, 
					{0, (2*pow(cw,2))/(pow(cw,2) - pow(sw,2)), 1/(pow(cw,2) - pow(sw,2)), 0, -1, 0, -((pow(cw,2)*pow(sw,2))/(pow(cw,2) - pow(sw,2))), 0, 0, 0}, 
					{0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, 
					{0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, 
					{0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, 
					{0, 0, 0, 0, 0, 0, 0, 0, 0, 1}};
 
    }
    else if (basis_output == bEFT_HiggsBasis){
      res.assign(nEFT_HiggsBasis_CouplingTypes, std::vector<double>(nEFT_HiggsBasis_CouplingTypes, 0));
      for (size_t i=0; i<(size_t) nEFT_HiggsBasis_CouplingTypes; i++) res.at(i).at(i)=1;
    }
  }

  if (res.empty()){
    cerr << "LexiConHCTranslator::getTranslationMatrix: Translation from input basis " << basis_input << " to output basis " << basis_output << " is not implemented." << endl;
    assert(0);
  }

  return res;
}

std::vector< std::pair<double, double> > LexiConHCTranslator::getOrderedInputCouplings(
  LexiConHCIOHelpers::IOBasisType const& basis_input,
  std::unordered_map<std::string, bool> const& input_flags,
  std::unordered_map<std::string, double> const& input_parameters,
  std::unordered_map<std::string, std::pair<double, double> > const& input_couplings
) const{
  bool useMCFMAtInput; getValueWithDefault<std::string, bool>(input_flags, "useMCFMAtInput", useMCFMAtInput, false);
  bool distinguish_HWWcouplings; getValueWithDefault<std::string, bool>(input_flags, "distinguish_HWWcouplings", distinguish_HWWcouplings, false);
  double MZ; getValueWithDefault<std::string, double>(input_parameters, "MZ", MZ, DEFVAL_MZ);
  double MW; getValueWithDefault<std::string, double>(input_parameters, "MW", MW, DEFVAL_MW);
  double Lambda_z1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_z1", Lambda_z1, DEFVAL_LAMBDA_VI);
  double Lambda_w1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_w1", Lambda_w1, DEFVAL_LAMBDA_VI);
  double Lambda_zgs1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_zgs1", Lambda_zgs1, DEFVAL_LAMBDA_VI);

  std::vector< std::pair<double, double> > res;
#define COUPLING_COMMAND(NAME, PREFIX, DEFVAL) \
  getValueWithDefault<std::string, std::pair<double, double>>(input_couplings, #NAME, res.at(coupl_##PREFIX##_##NAME), std::pair<double, double>(DEFVAL, 0)); \
  if (useMCFMAtInput && (std::string(#NAME).find("ghz")!=std::string::npos || std::string(#NAME).find("ghw")!=std::string::npos)){ res.at(coupl_##PREFIX##_##NAME).first *= 2.; res.at(coupl_##PREFIX##_##NAME).second *= 2.; } \
  if (std::string(#NAME).find("ghzgs")!=std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ res.at(coupl_##PREFIX##_##NAME).first *= std::pow(MZ/Lambda_zgs1, 2); res.at(coupl_##PREFIX##_##NAME).second *= std::pow(MZ/Lambda_zgs1, 2); } \
  else if (std::string(#NAME).find("ghz")!=std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ res.at(coupl_##PREFIX##_##NAME).first *= std::pow(MZ/Lambda_z1, 2); res.at(coupl_##PREFIX##_##NAME).second *= std::pow(MZ/Lambda_z1, 2); } \
  else if (std::string(#NAME).find("ghw")!=std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ res.at(coupl_##PREFIX##_##NAME).first *= std::pow(MW/Lambda_w1, 2); res.at(coupl_##PREFIX##_##NAME).second *= std::pow(MW/Lambda_w1, 2); }

  switch (basis_input){
  case bAmplitude_JHUGen:
  {
    res.assign(nAmplitude_JHUGen_CouplingTypes, std::pair<double, double>(0, 0));
    AMPLITUDE_JHUGEN_COUPLING_COMMANDS;
    break;
  }
  case bEFT_JHUGen:
  {
    res.assign(nEFT_JHUGen_CouplingTypes, std::pair<double, double>(0, 0));
    EFT_JHUGEN_COUPLING_COMMANDS;
    break;
  }
  case bHiggsBasis:
  {
    res.assign(nHiggsBasis_CouplingTypes, std::pair<double, double>(0, 0));
    HIGGSBASIS_COUPLING_COMMANDS;
    break;
  }
  case bEFT_HiggsBasis:
  {
    res.assign(nEFT_HiggsBasis_CouplingTypes, std::pair<double, double>(0, 0));
    EFT_HIGGSBASIS_COUPLING_COMMANDS;
    break;
  }
  default:
    cerr << "LexiConHCTranslator::getOrderedInputCouplings: Input basis " << basis_input << " is not implemented." << endl;
    assert(0);
  }

#undef COUPLING_COMMAND

  return res;
}

void LexiConHCTranslator::interpretOutputCouplings(
  LexiConHCIOHelpers::IOBasisType const& basis_output,
  std::unordered_map<std::string, bool> const& input_flags,
  std::unordered_map<std::string, double> const& input_parameters,
  std::vector< std::pair<double, double> >& output_vector
){
  bool useMCFMAtOutput; getValueWithDefault<std::string, bool>(input_flags, "useMCFMAtOutput", useMCFMAtOutput, false);

  double MZ; getValueWithDefault<std::string, double>(input_parameters, "MZ", MZ, DEFVAL_MZ);
  double MW; getValueWithDefault<std::string, double>(input_parameters, "MW", MW, DEFVAL_MW);
  double Lambda_z1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_z1", Lambda_z1, DEFVAL_LAMBDA_VI);
  double Lambda_w1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_w1", Lambda_w1, DEFVAL_LAMBDA_VI);
  double Lambda_zgs1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_zgs1", Lambda_zgs1, DEFVAL_LAMBDA_VI);
#define COUPLING_COMMAND(NAME, PREFIX, DEFVAL) \
  if (useMCFMAtOutput && (std::string(#NAME).find("ghz")!=std::string::npos || std::string(#NAME).find("ghw")!=std::string::npos)){ output_vector.at(coupl_##PREFIX##_##NAME).first /= 2.; output_vector.at(coupl_##PREFIX##_##NAME).second /= 2.; } \
  if (std::string(#NAME).find("ghzgs")!=std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ output_vector.at(coupl_##PREFIX##_##NAME).first /= std::pow(MZ/Lambda_zgs1, 2); output_vector.at(coupl_##PREFIX##_##NAME).second /= std::pow(MZ/Lambda_zgs1, 2); } \
  else if (std::string(#NAME).find("ghz")!=std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ output_vector.at(coupl_##PREFIX##_##NAME).first /= std::pow(MZ/Lambda_z1, 2); output_vector.at(coupl_##PREFIX##_##NAME).second /= std::pow(MZ/Lambda_z1, 2); } \
  else if (std::string(#NAME).find("ghw")!=std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ output_vector.at(coupl_##PREFIX##_##NAME).first /= std::pow(MW/Lambda_w1, 2); output_vector.at(coupl_##PREFIX##_##NAME).second /= std::pow(MW/Lambda_w1, 2); } \
  result_couplings[#NAME] = output_vector.at(coupl_##PREFIX##_##NAME);
  switch (basis_output){
  case bAmplitude_JHUGen:
  {
    AMPLITUDE_JHUGEN_COUPLING_COMMANDS;
    break;
  }
  case bEFT_JHUGen:
  {
    EFT_JHUGEN_COUPLING_COMMANDS;
    break;
  }
  case bHiggsBasis:
  {
    HIGGSBASIS_COUPLING_COMMANDS;
    break;
  }
  case bEFT_HiggsBasis:
  {
    EFT_HIGGSBASIS_COUPLING_COMMANDS;
    break;
  }
  default:
    cerr << "LexiConHCTranslator::interpretOutputCouplings: Output basis " << basis_output << " is not implemented." << endl;
    assert(0);
  }

#undef COUPLING_COMMAND
}
