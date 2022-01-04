/*
                           _ooOoo_
                          o8888888o
                          88" . "88
                          (| -_- |)
                          O\  =  /O
                       ____/`---'\____
                     .'  \\|     |//  `.
                    /  \\|||  :  |||//  \
                   /  _||||| -:- |||||-  \
                   |   | \\\  -  /// |   |
                   | \_|  ''\---/''  |   |
                   \  .-\__  `-`  ___/-. /
                 ___`. .'  /--.--\  `. . __
              ."" '<  `.___\_<|>_/___.'  >'"".
             | | :  `- \`.;`\ _ /`;.`/ - ` : | |
             \  \ `-.   \_ __\ /__ _/   .-` /  /
        ======`-.____`-.___\_____/___.-`____.-'======
                           `=---='
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                 佛祖保佑       永无BUG
        */



#include <iostream>
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include <stdio.h>
#include <fstream>
#include "TH1F.h"
#include "TClonesArray.h"
#include "coalescence.h"
#include  "TLegend.h"
#include <omp.h>
#define pi 3.14159
using namespace std;

int main(){

    
cout<<"                   _ooOoo_"<<endl;
cout<<"                 o8888888o"<<endl;
cout<<"      阿            88\" . \"88"<<endl;
cout<<"      弥            (| -_- |)"<<endl;
cout<<"      陀            O\\  =  /O"<<endl;
cout<<"               ____/\`---'\\____"<<endl;
cout<<"             .'  \\\\|  %  |//  `."<<endl;
cout<<"            /  \\\\|||  %  |||//  \\"<<endl;
cout<<"           /  _||||| -%- |||||-  \\"<<endl;
cout<<"           |   | \\\\\\  -  /// |   |"<<endl;
cout<<"           | \\_|  \'\'\\---/''  |   |"<<endl;
cout<<"           \\  .-\\__  `-`  ___/-. /"<<endl;
cout<<"         ___`. .'  /--.--\\  `. . __"<<endl;
cout<<"      .\"\" \'<  \`.___\\_<|>_/___.\'  >\'\"\"."<<endl;
cout<<"     | | :  \`- \\\`.;\`\\ _ /\`;.\`/ - \` : | |"<<endl;
cout<<"     \\  \\ \`-.   \\_ __\\ /__ _/   .-\` /  /"<<endl;
cout<<"======\`-.____\`-.___\\_____/___.-\`____.-\'======"<<endl;
cout<<"                   \`=---=\'"<<endl;
cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;
cout<<"              佛祖保佑        永无bug      "<<endl;








    const Int_t maxMultiplicity = 50000;
    Float_t b;
    Int_t nMultiplicityTree,npart;
    Int_t id[maxMultiplicity];
    Float_t px[maxMultiplicity],py[maxMultiplicity],pz[maxMultiplicity];
    Float_t x[maxMultiplicity],y[maxMultiplicity],z[maxMultiplicity],t[maxMultiplicity],mass[maxMultiplicity],energy[maxMultiplicity];

    TFile *hadronTreeFile = new TFile("onlyNeutron.root", "read");
    TTree *hadronTree = (TTree*)hadronTreeFile->Get("neutron");

	TH1F *teraneutron = new TH1F("teraneutron","teranuetron",40,0.,4);

    hadronTree->SetBranchAddress("mult",&nMultiplicityTree);
    hadronTree->SetBranchAddress("impact",&b);
    hadronTree->SetBranchAddress("pid",id);
    hadronTree->SetBranchAddress("fx",x);
    hadronTree->SetBranchAddress("fy",y);
    hadronTree->SetBranchAddress("fz",z);
    hadronTree->SetBranchAddress("fpx",px);
    hadronTree->SetBranchAddress("fpy",py);
    hadronTree->SetBranchAddress("fpz",pz);
    hadronTree->SetBranchAddress("ft",t);
    hadronTree->SetBranchAddress("fmass",mass);
    hadronTree->SetBranchAddress("fenergy",energy);

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

                if(coalescence(p1,p2,r1,r2,t[j],t[k])==false) continue;
                for(int l=k+1;l<nMultiplicityTree;l++){
                    //third particle
                    if(result==true)    break;
                    if(std::find(paired.begin(),paired.end(),l) !=paired.end())    continue;
                    TLorentzVector p3,r3;
                    p3.SetPxPyPzE(px[l],py[l],pz[l],energy[l]);
                    r3.SetPxPyPzE(x[l],y[l],z[l],t[l]);
                    if(coalescence(p1,p3,r1,r3,t[j],t[l])==false || coalescence(p2,p3,r2,r3,t[k],t[l])==false)  continue;
                    
                    #pragma omp parallel for
                    for(int m=l+1;m<nMultiplicityTree;m++){
                        //forth particle
                        if(result==true)    continue;
                        if(std::find(paired.begin(),paired.end(),m) !=paired.end())    continue;
                        TLorentzVector p4,r4;
                        p4.SetPxPyPzE(px[m],py[m],pz[m],energy[m]);
                        r4.SetPxPyPzE(x[m],y[m],z[m],t[m]);
                        if(coalescence(p1,p4,r1,r4,t[j],t[m])==true){
			
                        TLorentzVector four;
                        four=p1+p2+p3+p4;
                        float pt = four.Perp();
                        float rapidity = four.Rapidity();
                        if(rapidity>0.5 || rapidity<-0.5)   continue;
                        teraneutron->Fill(pt,1.0/(2.0*pi*0.1*2*pt));
                        paired.push_back(j);
                        paired.push_back(k);
                        paired.push_back(l);
                        paired.push_back(m);
                        result=true;
			continue;
			}
                        if(coalescence(p2,p4,r2,r4,t[k],t[m])==true){
                        TLorentzVector four;
                        four=p1+p2+p3+p4;
                        float pt = four.Perp();
                        float rapidity = four.Rapidity();
                        if(rapidity>0.5 || rapidity<-0.5)   continue;
                        teraneutron->Fill(pt,1.0/(2.0*pi*0.1*2*pt));
                        paired.push_back(j);
                        paired.push_back(k);
                        paired.push_back(l);
                        paired.push_back(m);
                        result=true;
			continue;
			}
                        if(coalescence(p3,p4,r3,r4,t[l],t[m])==true){
                        TLorentzVector four;
                        four=p1+p2+p3+p4;
                        float pt = four.Perp();
                        float rapidity = four.Rapidity();
                        if(rapidity>0.5 || rapidity<-0.5)   continue;
                        teraneutron->Fill(pt,1.0/(2.0*pi*0.1*pt));
                        paired.push_back(j);
                        paired.push_back(k);
                        paired.push_back(l);
                        paired.push_back(m);
                        result=true;
			continue;
			}

                    }
                }

            }
           
        }
	   paired.clear();
    }
    teraneutron->Scale(1.0/nentries);
    TCanvas *c1 = new TCanvas();
    teraneutron->Draw("e");
    teraneutron->GetXaxis()->SetTitle("p_{T} GeV");
    teraneutron->GetYaxis()->SetTitle("#frac{dN^{2}}{2#pip_{T}dp_{T}dy}");
    TLegend *leg1 = new TLegend();
    leg1->AddEntry(teraneutron,"Tetraneutron");
    c1->Draw();
    c1->SaveAs("tetraneutron.root");
    hadronTreeFile.Close();
    cout <<"done!" <<endl;    
    return 0;
}
