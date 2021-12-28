

void data2N(){

    const Int_t maxMultiplicity = 50000;
    Float_t b;
    Int_t nMultiplicityTree;
    Int_t id[maxMultiplicity];
    Float_t px[maxMultiplicity],py[maxMultiplicity],pz[maxMultiplicity];
    Float_t x[maxMultiplicity],y[maxMultiplicity],z[maxMultiplicity],t[maxMultiplicity],mass[maxMultiplicity],energy[maxMultiplicity];

    TFile *hadronTreeFile = new TFile("afterART2005.root", "read");
    TTree *hadronTree = (TTree*)hadronTreeFile->Get("hadronTree");

    hadronTree->SetBranchAddress("nMultiplicityTree",&nMultiplicityTree);
    hadronTree->SetBranchAddress("b",&b);
    hadronTree->SetBranchAddress("id",id);
    hadronTree->SetBranchAddress("x",x);
    hadronTree->SetBranchAddress("y",y);
    hadronTree->SetBranchAddress("z",z);
    hadronTree->SetBranchAddress("px",px);
    hadronTree->SetBranchAddress("py",py);
    hadronTree->SetBranchAddress("pz",pz);
    hadronTree->SetBranchAddress("t",t);
    hadronTree->SetBranchAddress("mass",mass);
    hadronTree->SetBranchAddress("energy",energy);

    //create a new tree to save data after decay.
    Int_t pid[maxMultiplicity];
    Int_t mult;
    Float_t impact;
    Float_t fpx[maxMultiplicity],fpy[maxMultiplicity],fpz[maxMultiplicity];
    Float_t fx[maxMultiplicity],fy[maxMultiplicity],fz[maxMultiplicity],ft[maxMultiplicity],fmass[maxMultiplicity],fenergy[maxMultiplicity];
    TFile *newFile = new TFile("onlyNeutron.root","RECREATE");
    TTree *neutron = new TTree("neutron","neutron"); 
    neutron->Branch("mult",&mult,"mult/I");
    neutron->Branch("impact",&impact,"impact/F");
    neutron->Branch("pid",&pid,"pid[mult]/I");
    neutron->Branch("fx",&fx,"fx[mult]/F");
    neutron->Branch("fy",&fy,"fy[mult]/F");
    neutron->Branch("fz",&fz,"fz[mult]/F");
    neutron->Branch("fpx",&fpx,"fpx[mult]/F");
    neutron->Branch("fpy",&fpy,"fpy[mult]/F");
    neutron->Branch("fpz",&fpz,"fpz[mult]/F");
    neutron->Branch("ft",&ft,"ft[mult]/F");
    neutron->Branch("fmass",&fmass,"fmass[mult]/F");
    neutron->Branch("fenergy",&fenergy,"fenergy[mult]/F");   


    const Int_t nevent = (Int_t)hadronTree->GetEntries();
    for(Int_t i=0;i<nevent;i++){
        hadronTree->GetEntry(nevent);
        Int_t m=0;
        impact=b;
        for(Int_t j=0;j<nMultiplicityTree;j++){
            if(id[j]!=2112)    continue; 
            pid[m]=id[j];
            fx[m]=x[j];
            fy[m]=y[j];
            fz[m]=z[j];
            fpx[m]=px[j];
            fpy[m]=py[j];
            fpz[m]=pz[j];
            ft[m]=t[j];
            fmass[m]=mass[j];
            fenergy[m]=energy[j];
            m++;
        }
        mult=m;
        neutron->Fill();
    }
    newFile->Write();
    cout<<"done!"<<endl;

}