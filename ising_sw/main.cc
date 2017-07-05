//----------------------------------------------------------------------
#include <iostream>
#include <math.h>
//----------------------------------------------------------------------
const int L = 64;       // System Size
const int T_LOOP = 100; // Thermalization loop
const int S_LOOP = 100; // Number of average
const int O_LOOP = 10;  // Number of observation at a temperature
const int ND=20;        // Number of temperatures to be observed 
int spin[L][L];
int cluster[L*L];
int flip[L*L];
//----------------------------------------------------------------------
inline double
myrand(void){
  return (double)rand()/(double)RAND_MAX;
}
//----------------------------------------------------------------------
void
init(void){
  for(int i=0;i<L;i++){
    for(int j=0;j<L;j++){
      spin[i][j] = 1;
    }
  }
}
//----------------------------------------------------------------------
int
get_cluster_number(int index){
  int i = index;
  while(i != cluster[i]){
    i = cluster[i];
  }
  return i;
}
//----------------------------------------------------------------------
void
connect_bond(int ix1, int iy1, int ix2, int iy2){
  const int i1 = ix1 + iy1*L;
  const int i2 = ix2 + iy2*L;
  const int c1 = get_cluster_number(i1);
  const int c2 = get_cluster_number(i2);
  if(c1<c2){
    cluster[c2] = c1;
  }else{
    cluster[c1] = c2;
  }
}
//----------------------------------------------------------------------
void
cluster_flip(double beta){
  const double p = 1.0 - exp(-2.0*beta);
  for(int i=0;i<L*L;i++){
    cluster[i] = i;
    flip[i] = static_cast<int>(myrand()*2.0) * 2 - 1;
  }

  for(int ix=0;ix<L-1;ix++){
    for(int iy=0;iy<L;iy++){
      if(spin[ix][iy] * spin[ix+1][iy] > 0 && myrand() < p){
        connect_bond(ix,iy,ix+1,iy);
      }
    }
  }

  for(int iy=0;iy<L;iy++){
    if(spin[L-1][iy] * spin[0][iy] > 0 && myrand() < p){
      connect_bond(L-1,iy,0,iy);
    } 
  } 

  for(int ix=0;ix<L;ix++){
    for(int iy=0;iy<L-1;iy++){
      if(spin[ix][iy] * spin[ix][iy+1] > 0 && myrand() < p){
        connect_bond(ix,iy,ix,iy+1);
      }
    }
  }

  for(int ix=0;ix<L;ix++){
    if(spin[ix][L-1] * spin[ix][0] > 0 && myrand() < p){
      connect_bond(ix,L-1,ix,0);
    } 
  } 

  for(int ix=0;ix<L;ix++){
    for(int iy=0;iy<L;iy++){
      const int c = get_cluster_number(ix+iy*L);
      spin[ix][iy] = spin[ix][iy] * flip[c];
    }
  }
}
//----------------------------------------------------------------------
void
single_flip(double beta){
  for(int ix=0;ix<L;ix++){
    for(int iy=0;iy<L;iy++){
      double de = 0;
      if(ix!=0)de += spin[ix-1][iy];
      else de += spin[L-1][iy];
      if(iy!=0)de += spin[ix][iy-1];
      else de += spin[ix][L-1];
      if(ix!=L-1)de += spin[ix+1][iy];
      else de += spin[0][iy];
      if(iy!=L-1)de += spin[ix][iy+1];
      else de += spin[ix][0];
      de = de*2.0*spin[ix][iy];
      if(de<0 || exp(-beta*de) > myrand()){
        spin[ix][iy] = - spin[ix][iy];
      }
    }
  }
}
//----------------------------------------------------------------------
void
mc_onestep(double beta){
#ifdef CLUSTER
  cluster_flip(beta);
#else
  single_flip(beta);
#endif
}
//----------------------------------------------------------------------
// Observation
//----------------------------------------------------------------------
double
magnetization(void){
  double m = 0;
  for(int ix=0;ix<L;ix++){
    for(int iy=0;iy<L;iy++){
      m += spin[ix][iy];
    }
  }
  m = m /(double)(L*L);
  if(m<0) m = -m;
  return m;
}
//----------------------------------------------------------------------
double
magnetization_square(void){
  double m = 0;
  for(int ix=0;ix<L;ix++){
    for(int iy=0;iy<L;iy++){
      m += spin[ix][iy];
    }
  }
  m /= static_cast<double>(L*L);
  return m*m;
}
//----------------------------------------------------------------------
double
energy(void){
  double e = 0.0;
  for(int ix = 0; ix< L;ix++){
    for(int iy = 0; iy< L;iy++){
      const int s = spin[ix][iy];
      if(ix!=L-1) e -= s*spin[ix+1][iy];
      else        e -= s*spin[0][iy];
      if(iy!=L-1) e -= s*spin[ix][iy+1];
      else        e -= s*spin[ix][0];
   }
  }
  e /= static_cast<double>(L*L);
  return e;
}
//----------------------------------------------------------------------
void
domc(void){
  double bs = 0.3; //Start inverse temperature  (The highest temperature)
  double be = 0.6; //End inverse temperature (The lowest temperature)
  for(int i=0;i<ND;i++){
    double beta = bs + (be-bs)*(double)i/(double)ND;
    fprintf(stderr,"beta = %f\n",beta);
    init();
    for(int j=0;j<T_LOOP;j++){
      mc_onestep(beta);
    }
    for(int j=0;j<O_LOOP;j++){
      double m2 = 0.0;
      double e_sum = 0.0;
      double e2_sum = 0.0;
      for(int k=0;k<S_LOOP;k++){
        mc_onestep(beta);
        m2 += magnetization_square();
        const double e = energy();
        e_sum += e;
        e2_sum += (e*e);
      }
      m2 /= static_cast<double>(S_LOOP);
      e_sum /= static_cast<double>(S_LOOP);
      e2_sum /= static_cast<double>(S_LOOP);
      const double c = (e2_sum - e_sum*e_sum)*beta*beta*L*L;
      printf("%f %f %f\n",beta,m2,c);
    }
  }
}
//----------------------------------------------------------------------
int
main(void){
#ifdef CLUSTER
  fprintf(stderr,"Swendsen-Wang\n");
  printf("# beta, Magnetization^2, Specific Heat\n");
#else
  fprintf(stderr,"Single Flip\n");
  printf("# beta, Magnetization^2, Specific Heat\n");
#endif
  domc();
}
//----------------------------------------------------------------------
