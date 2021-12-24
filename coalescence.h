
bool coalescence(TLorentzVector &p1,TLorentzVector &p2,TLorentzVector &r1,TLorentzVector &r2,float &t1,float &t2){
            
            TLorentzVector p3,p4,r3;              
            p3=p1+p2; //洛伦兹变换前的动量和
            if(fabs(p3.Rapidity())>0.3) return false;

            p1.Boost(-p3.BoostVector());
            p2.Boost(-p3.BoostVector());
            
            p4=p1-p2;
            if(p4.P()>0.22) return false;
          //  std::cout << "p4.p: "<<p4.P()<<std::endl;
            
            bool status = false;
    
            if(t1==t2)
            {
               
                    r1.Boost(-p3.BoostVector());
                    r2.Boost(-p3.BoostVector());        
                    r3=r1-r2;
                    if(r3.Rho()<2.3)   status = true;

                
            }
           else if(t1>t2) 
           {

                    Float_t dt=t1-t2;
                    TVector3 bv;
                    bv = p2.BoostVector();
                    TLorentzVector r22(bv.X()*dt+r2.X(),bv.Y()*dt+r2.Y(),bv.Z()*dt+r2.Z(),t1);
                    TLorentzVector deltaR;
                    deltaR=r1-r22;
                    deltaR.Boost(-p3.BoostVector());
                    if(deltaR.Rho()<2.3)    status = true;   
                    
           }
            
            else {
                    
                    Float_t dt=t2-t1;
                    TVector3 bv;
                    bv = p1.BoostVector();
                    TLorentzVector r11(bv.X()*dt+r1.X(),bv.Y()*dt+r1.Y(),bv.Z()*dt+r1.Z(),t2);                    
                    TLorentzVector deltaR;
                    deltaR=r11-r2;
                    deltaR.Boost(-p3.BoostVector());
                    if(deltaR.Rho()<2.3)    status = true; 

                }
      //      if(status==true)    std::cout << "Aha,paired!"<<std::endl;                       
        return status;
        
}