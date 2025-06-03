// system includes
#include <memory>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <functional>
#include <cassert>
#include <nlohmann/json.hpp>

//ROOT inclusion
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TMath.h"
#include "TList.h"
#include "TString.h"
#include "TLorentzVector.h"

using namespace std;

struct TreeReader{

  TTree* tree = NULL;
  
  long long nentries = 0;
  UInt_t event = 0;
  int run = 0;
  int lumi = 0;
  
  static constexpr int nMax = 40000;
  
  float hit_tob_lo_x[nMax] = {0},hit_tob_lo_y[nMax] = {0},hit_tob_lo_z[nMax] = {0};
  float hit_tib_lo_x[nMax] = {0},hit_tib_lo_y[nMax] = {0},hit_tib_lo_z[nMax] = {0};
  float hit_tid_lo_x[nMax] = {0},hit_tid_lo_y[nMax] = {0},hit_tid_lo_z[nMax] = {0};
  float hit_tec_lo_x[nMax] = {0},hit_tec_lo_y[nMax] = {0},hit_tec_lo_z[nMax] = {0};

  int hit_tob_detId[nMax] = {0},hit_tib_detId[nMax] = {0},hit_tid_detId[nMax] = {0},hit_tec_detId[nMax] = {0};

  float striphit_x=0,striphit_y=0,striphit_z=0;
  uint32_t detId = 0;
  
  TreeReader(TTree* in_tree):
    tree(in_tree)
  {
    nentries = tree->GetEntries();
    tree->SetBranchAddress("event",&event);
    tree->SetBranchAddress("run",&run);
    tree->SetBranchAddress("lumi",&lumi);
    tree->SetBranchAddress("detId",&detId);


    tree->SetBranchAddress("hit_tob_lo_x",hit_tob_lo_x);
    tree->SetBranchAddress("hit_tob_lo_y",hit_tob_lo_y);
    tree->SetBranchAddress("hit_tob_lo_z",hit_tob_lo_z);
    tree->SetBranchAddress("hit_tob_detId",hit_tob_detId);

    tree->SetBranchAddress("hit_tib_lo_x",hit_tib_lo_x);
    tree->SetBranchAddress("hit_tib_lo_y",hit_tib_lo_y);
    tree->SetBranchAddress("hit_tib_lo_z",hit_tib_lo_z);
    tree->SetBranchAddress("hit_tib_detId",hit_tib_detId);

    tree->SetBranchAddress("hit_tid_lo_x",hit_tid_lo_x);
    tree->SetBranchAddress("hit_tid_lo_y",hit_tid_lo_y);
    tree->SetBranchAddress("hit_tid_lo_z",hit_tid_lo_z);
    tree->SetBranchAddress("hit_tid_detId",hit_tid_detId);

    tree->SetBranchAddress("hit_tec_lo_x",hit_tec_lo_x);
    tree->SetBranchAddress("hit_tec_lo_y",hit_tec_lo_y);
    tree->SetBranchAddress("hit_tec_lo_z",hit_tec_lo_z);
    tree->SetBranchAddress("hit_tec_detId",hit_tec_detId);

    tree->SetBranchAddress("striphit_x",&striphit_x);
    tree->SetBranchAddress("striphit_y",&striphit_y);
    tree->SetBranchAddress("striphit_z",&striphit_z);

    
  };
  ~TreeReader(){
    delete tree;
  };
  
};

using json = nlohmann::json;

int main(int argc, char const *argv[]){
    
  TFile* f = TFile::Open(argv[1],"read");
  TreeReader tree( (TTree*)f->Get("clusterAnalyzer/clusterTree") );
  cout<<"The entries "<<tree.nentries<<endl;
  std::map <int,int> event_map;
  for (int ievt = 0; ievt < tree.nentries; ievt++){
    tree.tree->GetEntry(ievt);
    
    // int tobsize = sizeof(tree.hit_tob_detId);
    // for (int htob = 0; htob < tobsize; htob++){ 
    //   if(tree.hit_tob_detId[htob] == tree.detId){
    // 	TLorentzVector lz;
    // 	lz.SetXYZT(tree.hit_tob_lo_x[htob],tree.hit_tob_lo_y[htob],tree.hit_tob_lo_z[htob],0);
    //   }
    // }

    event_map[tree.event] = tree.detId;
  }
  // Convert std::map to nlohmann::json object
  json jsonObj = event_map;

  // Output JSON to console
  std::cout << jsonObj.dump(4) << std::endl;

  // Write JSON to a file
  std::ofstream outFile("output.json");

  if (outFile.is_open()) {
    outFile << jsonObj.dump(4); // Pretty-print with an indentation of 4 spaces
    outFile.close();
    std::cout << "JSON has been written to output.json" << std::endl;
  } else {
    std::cerr << "Failed to open file for writing." << std::endl;
  }

}
