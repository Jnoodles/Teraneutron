#include <iostream>
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include <stdio.h>
#include <fstream>
#include "TH1F.h"
#include "TClonesArray.h"
#include "coalescence.h"
#define pi 3.14159
using namespace std;

void tetraNeutron(){
    const Int_t maxMultiplicity = 50000;
    Float_t b;
    Int_t nMultiplicityTree,npart;
    Int_t id[maxMultiplicity];
    Float_t px[maxMultiplicity],py[maxMultiplicity],pz[maxMultiplicity];
    Float_t x[maxMultiplicity],y[maxMultiplicity],z[maxMultiplicity],t[maxMultiplicity],mass[maxMultiplicity],energy[maxMultiplicity];

    TFile *hadronTreeFile = new TFile("afterART2005.root", "read");
    TTree *hadronTree = (TTree*)hadronTreeFile->Get("hadronTree");

	TH1F *teraneutron = new TH1F("teraneutron","teranuetron",40,0.,2.0);

    hadronTree->SetBranchAddress("nMultiplicityTree",&nMultiplicityTree);
    hadronTree->SetBranchAddress("b",&b);
    hadronTree->SetBranchAddress("id",id);
    hadronTree->SetBranchAddress("x",x);
    hadronTree->SetBranchAddress("y",y);
    hadronTree->SetBranchAddress("z",z);
    hadronTree->SetBranchAddress("px",px);
    hadronTree->SetBranchAddress("py",py);
    hadronTree->SetBranchAddress("pz",pz);
    hadronTree->SetBranchAddress("t",&t);
    hadronTree->SetBranchAddress("mass",mass);
    hadronTree->SetBranchAddress("energy",energy);

	const Int_t nentries=hadronTree->GetEntries();

    cout<<"events: "<<nentries<<endl;

    vector<int> paired;

    for(int i=0;i<nentries;i++)
    {
		hadronTree->GetEntry(i);
    	for(int j=0;j<nMultiplicityTree;j++)
        {//first particle
            bool result =false;
            if(std::find(paired.begin(),paired.end(),j) !=paired.end())    continue;
            TLorentzVector p1,r1;
            p1.SetPxPyPzE(px[j],py[j],pz[j],energy[j]);
            r1.SetPxPyPzE(x[j],y[j],z[j],t[j]);

            for(int k=j+1;k<nMultiplicityTree;k++){
                //second particle
                if(result==true)    break;
                if(std::find(paired.begin(),paired.end(),k) !=paired.end())    continue;
                TLorentzVector p2,r2;
                p2.SetPxPyPzE(px[k],py[k],pz[k],energy[k]);
                r2.SetPxPyPzE(x[k],y[k],z[k],t[k]);

                if(coalescence(p1,p2,r1,r2)==false) continue;
                for(int l=k+1;l<nMultiplicityTree;l++){
                    //third particle
                    if(result==true)    break;
                    if(std::find(paired.begin(),paired.end(),l) !=paired.end())    continue;
                    TLorentzVector p3,r3;
                    p3.SetPxPyPzE(px[l],py[l],pz[l],energy[l]);
                    r3.SetPxPyPzE(x[l],y[l],z[l],t[l]);
                    if(coalescence(p1,p3,r1,r3)==false || coalescence(p2,p3,r2,r3)==false)  continue;
                    
                    for(int m=l+1;m<nMultiplicityTree;m++){
                        //forth particle
                        if(std::find(paired.begin(),paired.end(),m) !=paired.end())    continue;
                        TLorentzVector p4,r4;
                        p4.SetPxPyPzE(px[m],py[m],pz[m],energy[m]);
                        r4.SetPxPyPzE(x[m],y[m],z[m],t[m]);
                        if(colescence(p1,p4,r1,r4)==false) continue;
                        if(coalescence(p2,p4,r2,r4)==false) continue;
                        if(coalescence(p3,p4,r3,r4)==false) continue;

                        TLorentzVector four;
                        four=p1+p2+p3+p4;
                        float pt = four.Perp();
                        float rapidity = four.Rapidity();
                        if(rapidity>0.5 || rapidity<-0.5)   continue;
                        teraneutron->Fill(pt,1.0/(2.0*pi*0.05*pt));
                        paired.push_back(j);
                        paired.push_back(k);
                        paired.push_back(l);
                        paired.push_back(m);
                        result=true;
                        break;
                    }
                }

            }
           
        }
        paired.clear();
    }
    
    TCanvas *c1 = new TCanvas();
    teraneutron->Draw("e");
    c1->Draw();
    c1->SaveAs("tetraneutron.root");
    cout <<"done!" <<endl;    
}
