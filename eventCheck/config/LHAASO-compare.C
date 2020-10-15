{
#define CORE_X  0.0 //m    
#define CORE_Y  0.0 //m   
#define CORE_R  575.0 //m  
#define CORE_R2  637.0 //m 

  TCanvas *c1=new TCanvas("c1","The time offset effect",600,600);
  c1->SetFillColor(10);
  c1->SetFrameFillColor(kWhite);
 // c1->SetGridx(1);
//  c1->SetGridy(1);
//  c1->SetLogy(1);
//  c1->SetLogx(1);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.16);
  gStyle->SetStatStyle(0);
 gStyle->SetPalette(1); 
  gStyle->SetTitleXOffset(1.0);
  gStyle->SetTitleYOffset(1.);  

  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetOptStat(0);

  int i,j,id,k,nx,ny,n,ned,nmd;
  float dx,dy;
  TGraph *g[10];
  double x[6000],y[6000],z[6000],x2[6000],y2[6000],z2[6000];
   double xx,yy,zz;
   char sTemp[200];
   n=0;
     x[0]=-700;y[0]=-700;
     x[1]=700;y[1]=700;
     g[0]=new TGraph(2,x,y);

  for(i=0;i<6000;i++){
     x[i]=-1000; y[i]=0;z[i]=0;
     x2[i]=-1000; y2[i]=0;z2[i]=0;
  }
  ifstream file1,file2;
      file1.open("ED_pos_quarter.txt");
      file1.getline(sTemp,200);
     while(file1.good()){
           file1.getline(sTemp,200);
           i=sscanf(sTemp,"%d %lf %lf %lf",&id,&xx,&yy,&zz);
           if(i>1){
              n++;
              x[id]=xx; y[id]=yy;z[id]=zz;
           }
      }
      file1.close();
   ned=n;
   g[1]=new TGraph(6000,x,y);
   printf("No. of ED is %d\n",n);

   n=0;
   //file1.open("MD_pos_half.txt");
   file1.open("/scratchfs/ybj/zhanghy/KM2A/G4KM2A4.10/config/ED_pos_quarter.txt");
   // file1.open("../../test/G4KM2Arec6/config/ED_pos_quarter.txt");
   file1.getline(sTemp,200);
      while(file1.good()){
           file1.getline(sTemp,200);
           i=sscanf(sTemp,"%d %lf %lf %lf",&id,&xx,&yy,&zz);
           if(i>1){
              n++;
              x2[id]=xx; y2[id]=yy;z2[id]=zz;
           }
      }
      file1.close();
   nmd=n;
    g[2]=new TGraph(6000,x2,y2);
   printf("No. of ED is %d\n",n);

     for(i=0;i<6000;i++){
        if(x[i]>-990){
          if(sqrt(pow(x[i]-x2[i],2)+pow(y[i]-y2[i],2)+pow(z[i]-z2[i],2))>1){
             printf("%d %lf %lf %lf\n",i,x[i]-x2[i],y[i]-y2[i],z[i]-z2[i]);
         
          }
        }
     }

     g[0]->Draw("AP");
    if(ned>0) g[1]->Draw("P"); 
    if(nmd>0) g[2]->Draw("P");
    for(i=1;i<3;i++){
       g[i]->SetMarkerStyle(22-i);
       g[i]->SetMarkerColor(i+1);
       g[i]->SetMarkerSize(0.3);
    }
   g[1]->SetMarkerColor(2);
   g[2]->SetMarkerColor(4);
    g[0]->SetTitle(";X (m);Y (m)");
     g[0]->GetYaxis()->CenterTitle(true);
     g[0]->GetYaxis()->SetTitleOffset(1.1);
     g[0]->GetXaxis()->SetTitleOffset(1.1);
     g[0]->GetYaxis()->SetLabelSize(0.03);
     g[0]->GetXaxis()->SetLabelSize(0.03);
     g[0]->GetXaxis()->CenterTitle(true);
     g[0]->GetXaxis()->SetRangeUser(-650,650);
     g[0]->GetYaxis()->SetRangeUser(-650,650); 
     c1->Update();
     sprintf(sTemp,"LHAASO Layout: %d EDs + %d MDs",ned,nmd);
        TLatex *   tex = new TLatex(-500,700,sTemp);
        tex->SetLineWidth(2);
        tex->SetTextSize(0.04);
        tex->SetTextColor(6);
        tex->Draw();

      TEllipse *ellipse = new TEllipse(CORE_X,CORE_Y,CORE_R+10,CORE_R+10,0,360,0);ellipse->SetFillStyle(0);
      // ellipse->SetLineWidth(2);
      ellipse->SetLineColor(13);
      ellipse->Draw();  

      ellipse = new TEllipse(CORE_X,CORE_Y,CORE_R2,CORE_R2,0,360,0);ellipse->SetFillStyle(0);
      //ellipse->SetLineWidth(2);
      ellipse->SetLineColor(13);
      ellipse->Draw();
 
    c1->Update();
  return 0;
}
